import os
import json
import re
import subprocess
import xmltodict
import numpy as np
from Bio.PDB.DSSP import DSSP
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import MMCIFIO, MMCIFParser
from Bio.PDB import NeighborSearch
from Bio import SwissProt
from Bio.PDB.Polypeptide import PPBuilder
import freesasa
import networkx as nx
import getBASA
from dnaprodb_utils import getID, getHash, log, C
from dnaprodb_utils import CHAIN_RE, RESN_RE, RESI_RE, ATOM_RE
from dnaprodb_utils import getStructureFromModel

# External data directories
DATA_PATH = C["DATA_PATH"]
ROOT_DIR = C["ROOT_DIR"]
UNIPROT_DIR = os.path.join(ROOT_DIR, "UNIPROT")
PDB_DIR = os.path.join(ROOT_DIR, "PDB")
CATH_DIR = os.path.join(ROOT_DIR, "CATH")
SIFTS_DIR = os.path.join(ROOT_DIR, "SIFTS")

# Import standard SASA values
with open(os.path.join(DATA_PATH,'standard-sasa.json')) as FILE:
    STANDARD_SASA = json.load(FILE)

# Import residue hydrophobicities
with open(os.path.join(DATA_PATH,'residue-hydrophobicity.json')) as FILE:
    RESIDUE_HPHOB = json.load(FILE)

# Get various cut-off values
SAP_CUTOFF_DISTANCE = C["SAP_CUTOFF_DISTANCE"]
CV_FINE_CUTOFF_DISTANCE = C["CV_FINE_CUTOFF_DISTANCE"]
CV_COARSE_CUTOFF_DISTANCE = C["CV_COARSE_CUTOFF_DISTANCE"]

# Load additional data
shortNameMap = C["LONG_TO_SHORT"]
DSSP_MAP = C["DSSP_MAP"]

def getSESA(model, residue_ids, classifier, pdbid, regexes=None, total=True):
    FH = open("coords.xyzr","w")
    resData = []
    residueSESA = {} # keyed by residue ID
    for rid in residue_ids:
        cid, resi, ins = rid.split('.')
        residue = model[cid][(' ', int(resi), ins)]
        resn = residue.get_resname()
        residueSESA[rid] = {
            'sc': 0.0,
            'mc': 0.0,
            'total': 0.0
        }
        for atom in residue.get_list(): 
            if(atom.element == 'H'):
                continue
            aname = atom.name
            acoords = atom.get_coord()
            radius = classifier.radius(resn, aname)
            FH.write("{}\t{}\t{}\t{}\n".format(acoords[0], acoords[1], acoords[2],radius))
            resData.append((rid,resn,aname))
    FH.close()
    
    # Call mslib and compute sesa
    args = ["msms", "-surface", "ases", "-probe_radius", "1.4", "-if", "coords.xyzr", "-af", "mol.area"]
    FNULL = open(os.devnull, 'w')
    subprocess.call(args, stdout=FNULL, stderr=FNULL)
    FNULL.close()
    
    # Read in atomic SESA
    SE = open("mol.area").readlines()[1:]
    os.remove("mol.area")
    os.remove("coords.xyzr")
    if(len(SE) != len(resData)):
        log("Area and coordinate file lengths do not match in getSESA()!", pdbid)
    
    if(total):
        totalSESA = 0.0
        for i in range(len(SE)):
            totalSESA += float(SE[i].strip().split()[1])
        return totalSESA
    else:
        for i in range(len(SE)):
            rid, resn, aname = resData[i]
            sesa = float(SE[i].strip().split()[1])
            residueSESA[rid]["total"] += sesa
            if(regexes[resn][0].search(aname)):
                residueSESA[rid]["mc"] += sesa
            else:
                residueSESA[rid]["sc"] += sesa
        
        return residueSESA

def computeCV(model, res_list, residues, radius):
    # list of atom objects in the model
    atoms = []
    for rid in res_list:
        chain, num, ins = rid.split(".")
        residue = model[chain][(' ', int(num), ins)]
        for atom in residue.get_list(): 
            if(atom.element == 'H'):
                continue
            atoms.append(atom)
    ns = NeighborSearch(atoms)
    
    # get CV values for every residue in residues
    residueCV = []
    for resID in residues:
        chain, num, ins = resID.split('.')
        residue = model[chain][(' ', int(num), ins)]
        cvs = []
        for atom in residue:
            if(atom.element == 'H'):
                continue
            vector = np.zeros(3)
            neighbors = ns.search(atom.get_coord(), radius, level='A')
            for n in neighbors:
                if(n == atom):
                    continue
                vector += (atom.get_coord() - n.get_coord())/np.linalg.norm(atom.get_coord() - n.get_coord())
            cvs.append(1 - np.linalg.norm(vector)/(len(neighbors)-1))
        residueCV.append({"id": resID, "cv": np.mean(cvs)})
    
    return residueCV

def runDSSP(model, fileName):
    '''
    print(fileName)
    io=MMCIFIO()
    io.set_structure(model)
    cifname = fileName.replace(".pdb","_dssp.cif")
    io.save(cifname)
    
    cifparser = MMCIFParser()
    model = cifparser.get_structure(cifname, cifname)
    print(model, cifname)
    '''
    dssp = DSSP(model, fileName)
    ssMap = {}
    for chain in model:
        cid = chain.get_id()
        for residue in chain.get_list():
            rid = residue.get_id()
            dkey = (cid, rid)
            skey = getID(residue=residue)
            if(dkey in dssp):
                info = dssp[dkey]
                ssMap[skey] = DSSP_MAP[info[2]]
            else:
                ssMap[skey] = 'L'
    return ssMap

def getDSSP(arg, structure, fileName, chain_map=None):
    if(arg == 'consensus'):
        dssp = {}
        SS = C["RES_SST"]
        for model in structure:
            m = runDSSP(model, fileName)
            for key in m:
                if key in dssp:
                    dssp[key].append(m[key])
                else:
                    dssp[key] = [m[key]]
        
        for key in dssp:
            imax = 0
            cmax = 0
            for i in range(len(SS)):
                c = dssp[key].count(SS[i])
                if(c >= cmax):
                    imax = i
                    cmax = c
            dssp[key] = SS[imax]
        return dssp
    elif(type(arg) is int or arg.isdigit()):
        model = structure[int(arg)]
        return runDSSP(model, fileName)
    else:
        dssp = DSSP(structure[0], fileName)
        ssMap = {}
        
        for chain in structure[0]:
            cid = chain.get_id()
            if(chain_map):
                if(cid in chain_map):
                    ch = chain_map[cid]
                else:
                    ch = cid
            else:
                ch = cid
            for residue in chain.get_list():
                rid = residue.get_id()
                dkey = (cid, rid)
                skey = "{}.{}.{}".format(ch, rid[1], rid[2])
                if(dkey in dssp):
                    info = dssp[dkey]
                    ssMap[skey] = DSSP_MAP[info[2]]
                else:
                    ssMap[skey] = 'L'
        return ssMap

def addExternalData(models, pdbid, mmcif_dict, chain_map):
    clusters = ["30", "40", "50", "70", "90", "95", "100"]
    prefix = pdbid[-1]
    
    # Get UNIPROT data
    UNP_ACCS = {}
    UNP_RECORD = {}
    UNP_mapping = os.path.join(SIFTS_DIR, prefix, 'uniprot_mapping.tsv')
    if(os.access(UNP_mapping, os.R_OK)):
        FH = open(UNP_mapping)
        for line in FH:
            line = line.split()
            if(line[0] == pdbid):
                UNP_ACCS[line[1]] = line[2]
    
    for i in range(len(mmcif_dict['_struct_ref_seq.pdbx_strand_id'])):
        accession = None
        UNP_file = None
        if(mmcif_dict['_struct_ref_seq.pdbx_strand_id'][i] in UNP_ACCS):
            accession = UNP_ACCS[mmcif_dict['_struct_ref_seq.pdbx_strand_id'][i]]
            UNP_file = os.path.join(UNIPROT_DIR, accession[-1], accession+".txt")
        
        if(accession is None or not os.access(UNP_file, os.R_OK)):
            # Try another accession
            accession = mmcif_dict['_struct_ref_seq.pdbx_db_accession'][i]
            UNP_file = os.path.join(UNIPROT_DIR, accession[-1], accession+".txt")
        
        if(accession is None or not os.access(UNP_file, os.R_OK)):
            # Can't find an existing accession for this chain
            continue
        
        if(accession not in UNP_RECORD):
            handle = open(UNP_file)
            UNP_RECORD[accession] = SwissProt.read(handle)
            handle.close()
            names = []
            description = UNP_RECORD[accession].description
            description = re.sub(r'{.*?}', '', description)
            description = re.split(':|;',description)
            for i in range(len(description)):
                description[i] = description[i].strip()
                if(re.search('^Full|^Short',description[i])):
                    names.append(description[i].split('=')[1])
            UNP_RECORD[accession].description = names
    
    # Get CATH Data
    search = False
    if(os.access(os.path.join(CATH_DIR,'cath-domain-list.txt'), os.R_OK)):
        try:
            search = subprocess.check_output(['grep', pdbid, os.path.join(CATH_DIR,'cath-domain-list.txt')])
        except subprocess.CalledProcessError as e:
            pass
    CATH_DOMAINS = {}
    if(search):
        search = search.strip('\n').split('\n')
        for line in search:
            line = line.split()
            chain = line[0][4]
            cath_H = '.'.join(line[1:5])
            cath_T = '.'.join(line[1:4])
            cath_A = '.'.join(line[1:3])
            cath_C = line[1]
            if(chain in CATH_DOMAINS):
                if(cath_H not in CATH_DOMAINS[chain]['H']):
                    CATH_DOMAINS[chain]['H'].append(cath_H)
                if(cath_T not in CATH_DOMAINS[chain]['T']):
                    CATH_DOMAINS[chain]['T'].append(cath_T)
                if(cath_A not in CATH_DOMAINS[chain]['A']):
                    CATH_DOMAINS[chain]['A'].append(cath_A)
                if(cath_C not in CATH_DOMAINS[chain]['C']):
                    CATH_DOMAINS[chain]['C'].append(cath_C)
            else:
                CATH_DOMAINS[chain] = {
                    'H': [cath_H],
                    'T': [cath_T],
                    'A': [cath_A],
                    'C': [cath_C]
                }
    
    # Add data to each chain
    for model in models:
        for chain in model["chains"]:
            cid = chain['au_chain_id']
            for cluster in clusters:
                # Get sequence clusters from PDB
                if(os.access("{}/{}/{}.{}_{}.xml".format(PDB_DIR, prefix, pdbid, cid, cluster), os.R_OK)):
                    REP = open("{}/{}/{}.{}_{}.xml".format(PDB_DIR, prefix, pdbid, cid, cluster))
                    data = xmltodict.parse(REP.read())
                    REP.close()
                    if(data['representatives']):
                        chain['sequence_clusters'][cluster] = data['representatives']['pdbChain']['@name']
                    else:
                        chain['sequence_clusters'][cluster] = 'N/A'
                else:
                    chain['sequence_clusters'][cluster] = 'N/A'
            # Add UNIPROT
            if(cid in UNP_ACCS):
                accession = UNP_ACCS[cid]
            elif(chain_map[cid] in UNP_ACCS):
                accession = UNP_ACCS[chain_map[cid]]
            else:
                accession = None
            
            if(accession in UNP_RECORD):
                record = UNP_RECORD[accession]
                chain['uniprot_accession'] = record.accessions
                chain['uniprot_names'] = record.description
                chain['organism'] = record.organism
                chain["GO_molecular_function"] = []
                chain["GO_biological_process"] = []
                chain["GO_cellular_component"] = []
                for DR in record.cross_references:
                    if(DR[0] == 'GO'):
                        if(DR[2][0] == "F"):
                            chain["GO_molecular_function"].append({
                                "GO_ID": DR[1][3:],
                                "description": DR[2][2:]
                            })
                        elif(DR[2][0] == "P"):
                            chain["GO_biological_process"].append({
                                "GO_ID": DR[1][3:],
                                "description": DR[2][2:]
                            })
                        elif(DR[2][0] == "C"):
                            chain["GO_cellular_component"].append({
                                "GO_ID": DR[1][3:],
                                "description": DR[2][2:]
                            })
            else:
                chain['uniprot_accession'] = ['N/A']
                chain['uniprot_names'] = ['N/A']
                chain['organism'] = 'N/A'
                chain["GO_molecular_function"] = [{"description": 'N/A', "GO_ID": 'N/A'}]
                chain["GO_biological_process"] = [{"description": 'N/A', "GO_ID": 'N/A'}]
                chain["GO_cellular_component"] = [{"description": 'N/A', "GO_ID": 'N/A'}]
            # Add CATH
            if(cid in CATH_DOMAINS):
                chain['cath_homologous_superfamily'] = CATH_DOMAINS[cid]['H']
                chain['cath_topology'] = CATH_DOMAINS[cid]['T']
                chain['cath_architecture'] = CATH_DOMAINS[cid]['A']
                chain['cath_class'] = CATH_DOMAINS[cid]['C']
            else:
                chain['cath_homologous_superfamily'] = ['N/A']
                chain['cath_topology'] = ['N/A']
                chain['cath_architecture'] = ['N/A']
                chain['cath_class'] = ['N/A']

def _addLoops(index, sse, peptide, Lcount, SS_ELEMENTS):
    for j in range(sse["start"], index):
        # Add loops
        Lcount += 1
        loopID = peptide['residue_ids'][j]
        SS_ELEMENTS[loopID] = {
            "id": loopID,
            "secondary_structure": "L",
            "chain": peptide["chain"],
            "residue_ids": [peptide['residue_ids'][j]],
            "sequence": peptide['sequence'][j],
            "number": Lcount
        }
        peptide["secondary_structure"][j] = "L"
    return Lcount

def _addHelixStrand(index, sse, peptide, Count, SS_ELEMENTS):
    Count += 1
    sseID = '{}{}{}'.format(sse["type"], peptide["chain"], Count)
    SS_ELEMENTS[sseID] = {
        "id": sseID, 
        "secondary_structure": sse["type"], 
        "chain": peptide["chain"],
        "residue_ids": peptide['residue_ids'][sse["start"]:index],
        "sequence": peptide['sequence'][sse["start"]:index],
        "number": Count
    }
    return Count, sseID

def getPeptides(chain_info, chain):
    chain_info["id"]
    PPB = PPBuilder()
    peptides = PPB.build_peptides(chain, aa_only=False)
    PEPTIDES = []
    count = 1
    for pep in peptides:
        PEP = {
            "sequence": "",
            "secondary_structure": "",
            "residue_ids": [],
            "length": 0,
            "chain": chain_info["id"],
            "id": chain_info["id"]+str(count),
        }
        for r in pep:
            rid = getID(residue=r)
            index = chain_info["residue_ids"].index(rid)
            PEP["residue_ids"].append(rid)
            PEP["sequence"] += chain_info["sequence"][index]
            PEP["secondary_structure"] += chain_info["secondary_structure"][index]
            PEP["length"] += 1
        PEPTIDES.append(PEP)
        count += 1
    if(len(PEPTIDES) > 1):
        chain_info["continuous"] = False
    return PEPTIDES

def getSSE(peptides, structure, chains, residues):
    Hnum = 0
    Snum = 0
    HELIX = []
    SHEET = []
    SS_ELEMENTS = {}
    HcountCH = {}
    ScountCH = {}
    LcountCH = {}
    for peptide in peptides:
        ss = peptide['secondary_structure'] + "X"
        peptide["secondary_structure"] = list(peptide["secondary_structure"])
        cid = peptide["chain"]
        Hcount = HcountCH.get(cid, 0)
        Scount = ScountCH.get(cid, 0)
        Lcount = LcountCH.get(cid, 0)
        pre = ss[0]
        sse = {
            "start": 0,
            "type": ss[0]
        }
        for i in range(1, len(ss)):
            cur = ss[i]
            if(cur != pre):
                if(sse["type"] == "L"):
                    Lcount = _addLoops(i, sse, peptide, Lcount, SS_ELEMENTS)
                else:
                    length = i - sse["start"]
                    if(length >= 2):
                        # Add helix/strand
                        if(sse["type"] == "H"):
                            Hcount, helixID = _addHelixStrand(i, sse, peptide, Hcount, SS_ELEMENTS)
                            Hnum += 1
                            HELIX.append(printSSRemark(peptide['residue_ids'][sse["start"]:i], Hnum, helixID, structure,'H'))
                        else:
                            Scount, strandID = _addHelixStrand(i, sse, peptide, Scount, SS_ELEMENTS)
                            Snum += 1
                            SHEET.append(printSSRemark(peptide['residue_ids'][sse["start"]:i], Snum, strandID, structure,'S'))
                    else:
                        Lcount = _addLoops(i, sse, peptide, Lcount, SS_ELEMENTS)
                # make a new sse
                sse["start"] = i
                sse["type"]  = cur
            pre = cur
        peptide["secondary_structure"] = "".join(peptide["secondary_structure"])
        HcountCH[cid] = Hcount
        ScountCH[cid] = Scount
        LcountCH[cid] = Lcount
    
    # Ensure consistency between SSE and chain/residue secondary structure
    for chain in chains:
        chain["secondary_structure"] = list(chain["secondary_structure"])
        for i in range(len(chain["residue_ids"])):
            rid = chain["residue_ids"][i]
            if(rid in SS_ELEMENTS):
                chain["secondary_structure"][i] = "L"
                residues[rid]["secondary_structure"] = "L"
        chain["secondary_structure"] = "".join(chain["secondary_structure"])
    
    return list(SS_ELEMENTS.values()), HELIX, SHEET

def printSSRemark(resIDs, num, label, structure, sst):
    startID = resIDs[0]
    startChain, startNum, startIns = startID.split('.')
    startName = structure[0][startChain][(' ', int(startNum), startIns)].get_resname()
    
    endID = resIDs[-1]
    endChain, endNum, endIns = endID.split('.')
    endName = structure[0][endChain][(' ', int(endNum), endIns)].get_resname()
    
    if(sst == 'H'):
        return "{:<6s} {:>3d} {:>3s} {:>3s} {} {:>4s}{} {:>3s} {} {:>4s}{}{:>2d}{:>30s} {:>5d}".format(
            'HELIX',num,label[1:],startName,startChain,startNum,startIns,
            endName,endChain,endNum,endIns,1,"",len(resIDs))
    elif(sst == 'S'):
        return "{:<6s} {:>3d} {:>3s}{:>2d} {:>3s} {}{:>4s}{} {:>3s} {}{:>4s}{}{:>2d}{:>30s}".format(
            'SHEET',num,label[1:],0,startName,startChain,startNum,startIns,
            endName,endChain,endNum,endIns,0,"")

def stripH(model):
    """Strip all hydrogen atoms from the given model.
    
    Parameters
    ----------
    model: BioPython MODEL object
        The model to be stripped of hydrogens.
    """
    for chain in model:
        for residue in chain:
            rm = []
            for atom in residue:
                if(atom.element == 'H'):
                    rm.append(atom.get_id())
            for aid in rm:
                residue.detach_child(aid)

def getSurfaceResidues(res_ids, model=None, REGEXES=None, IDs=None, method='sesa', sesa=None):
    surface_residues= []
    
    if(method == 'sasa'):
        STRUCT = getBASA.buildStructure(model, res_ids)
        SASA = freesasa.calc(STRUCT)
        RES_SASA = getBASA.sumSASA(STRUCT, SASA, REGEXES, res_ids)
        for resID in RES_SASA:
            if(RES_SASA[resID]["sasa"][0]/STANDARD_SASA[RES_SASA[resID]["resn"]]["total"] > 0.05):
                surface_residues.append(resID)
    else:
        for resID in sesa:
            if(sesa[resID]["total"] >= 1.0):
                surface_residues.append(resID)
    
    return surface_residues

def getSAP(model, distance, classifier, REGEXES):
    # Reads in the file pdbid-protein.pdb and computes the SAP score for
    # all standard residues on the protein surface. Non-standard should 
    # be assigned a score of none. Non-standard residue atoms are ignored
    # and are not considered in SAP calculations.
    # Arguments:
    # pdbid:       structure name/identifier
    # distance:    distance cut-off for neighbor search
    #-------------------------------------------------------------------

    F_pro = getStructureFromModel(model, classifier=classifier)
    SASA = freesasa.calc(F_pro)
    
    # store all atoms with a SASA > 0, skipping non-standard atoms
    atomSASA = {}
    N = F_pro.nAtoms()
    for i in range(N):
        sasa = SASA.atomArea(i)
        if(sasa > 0 and F_pro.residueName(i).strip() in RESIDUE_HPHOB):
            resn = F_pro.residueName(i).strip()
            resi = F_pro.residueNumber(i).strip()
            if(resi[-1].isdigit()):
                ins = " "
            else:
                ins = resi[-1]
                resi = resi[:-1]
            chain = F_pro.chainLabel(i)
            rid = getID(chain, resi, ins)
            aname = F_pro.atomName(i).strip()
            atm_id = getID(rid, aname)
            atomSASA[atm_id] = sasa
    
    # list of atom objects corresponding to the atoms in atomSASA
    atoms = []
    for chain in model.get_list():
        cid = chain.get_id()
        for residue in chain.get_list():
            resn = residue.get_resname().strip()
            rid = getID(residue=residue)
            for atom in residue.get_list():
                aname = atom.get_name().strip()
                atm_id = getID(rid, aname) 
                if(atm_id in atomSASA):
                    atoms.append(atom)
    ns = NeighborSearch(atoms)
    
    # compute the residue SAP scores
    resSAP = {}
    for a in atoms:
        center = a.get_coord()
        res = a.get_parent()
        res_id = getID(residue=res)
        resn = res.get_resname().strip()
        rid = res.get_id()
        chain = res.get_parent()
        cid = chain.get_id()
        asap = 0
        neighbors = ns.search(center, distance, level='A')
        for n in neighbors:
            nname = n.get_name().strip()
            nres = n.get_parent()
            nrid = getID(residue=nres)
            nresn = nres.get_resname().strip()
            n_id = getID(nrid,nname)
            if(REGEXES['PROTEIN']['SIDE_CHAIN'].search(nname)):
                # select side-chain atoms only
                asap += RESIDUE_HPHOB[nresn]*(atomSASA[n_id])/STANDARD_SASA[nresn][nname]
        if(res_id not in resSAP):
            resSAP[res_id] = {
                'resn': resn,
                'resi': rid[1],
                'ins': rid[2],
                'chain': cid,
                'id': res_id,
                'acount': 0.0,
                'SAP': 0.0
            }
        resSAP[res_id]['acount'] += 1
        resSAP[res_id]['SAP'] += asap
    
    for key in resSAP:
        resSAP[key]['SAP'] /= resSAP[key]['acount']
    
    return list(resSAP.values())

def process(prefix, N, COMPONENTS, REGEXES, IDs,
    quiet=True,
    mmcif_dict=None,
    sse_type='per_model',
    reference_structure=None,
    chain_map=None,
    meta_data=None
    ):
    """Main function which processes a protein structure and extracts
    various information from it. Each model in the structure is stored
    as a JSON object in an array.
    
    Parameters
    ----------
    prefix: string
        Filename prefix string of a PDB file which contains only 
        protein residues and has all ligands, DNA, and solvent removed.
    N: int
        Number of models to process.
    quiet: boolean (True)
        Set to true to suppress output from external programs.
    annotate: boolean (False)
        If true, then annotate each chain with data from external 
        databases.
    sse_type: string | int ('per model')
        Decide how secondary structure is determined. Possible options 
        are the following:
        per_model - compute secondary structure separately for each 
            model in the structure. 
        consensus - compute secondary structure based on the consenus 
            across all models in the structure. This is redundant if 
            only one model is present.
        reference - assign secondary structure based on a reference 
            structure. This reference structure should correspond 
            exactly to the target structure except for the atomic 
            coordinates. If chosen, a reference structure must be 
            provided or an error will be generated.
        integer - if an integer is provided, a specific model number 
            will be used to assign secondary structure to all other 
            models.
    reference_structure: string (None)
        Needed if 'reference' is chosen for sse_type. This is the file 
        name of a reference PDB or mmCIF file.
    chain_map: dict (None)
        A mapping between the chains of the reference_structure and the 
        chains of the target structure. Only used for sse_type='reference'.
    meta_data: dict (None)
        The meta data dictionary passed from processStructure which contains
        various info about the structure processing.
    """
    
    fileName = "{}-protein.pdb".format(prefix)
    parser = PDBParser(PERMISSIVE=1, QUIET=quiet)
    structure = parser.get_structure("protein", fileName)
    
    # Add meta data if not provided
    if(meta_data is None):
        # Fill in missing data if meta_file not found
        meta_data = {
            "assembly_chains": {}
        }
        for chain in structure[0]:
            cid = chain.get_id()
            meta_data["assembly_chains"][cid] = cid
    
    # Create freeSASA classifier
    if(quiet):
        freesasa.setVerbosity(freesasa.nowarnings)
    classifier = getBASA.GeneralClassifier()
    classifier.initialize(DATA_PATH, COMPONENTS)
    
    # Output array - one entry for each model
    OUT = []
    MODEL_DSSP = [] # store dssp map for each model
    if(sse_type == 'reference'):
        reference = parser.get_structure("reference", reference_structure)
        dssp = getDSSP(sse_type, reference, reference_structure, chain_map=chain_map)
    elif(sse_type != 'per_model'):
        dssp = getDSSP(sse_type, structure, fileName)
    
    for i in range(N):
        # Get Data for each chain
        PRO_CHAINS = []
        PRO_RES = []
        model = structure[i]
        if(sse_type == 'per_model'):
            dssp = runDSSP(model, fileName)
        MODEL_DSSP.append(dssp)
        for chain in model:
            cid = chain.get_id()
            sst = "" # secondary structure sequence
            seq = "" # primary sequence
            ids = []
            for residue in chain.get_list():
                rname = residue.get_resname().strip()
                rid = residue.get_id()
                if(rname in shortNameMap):
                    seq += shortNameMap[rname]
                else:
                    seq += 'X'
                dkey = "{}.{}.{}".format(cid, rid[1], rid[2])
                if(dkey in dssp):
                    ss = dssp[dkey]
                else:
                    ss = 'L'
                sst += ss
                ids.append(dkey)
                PRO_RES.append({
                    "name": rname,
                    "name_short": shortNameMap.get(rname, "X"),
                    "chain": cid,
                    "secondary_structure": ss,
                    "id": dkey,
                    "ins_code": rid[2],
                    "number": rid[1],
                    "fasa":None,
                    "sesa": None,
                    "sap_score": None,
                    "cv_fine": None,
                    "cv_coarse": None,
                    "surface": False,
                    "chemical_name": COMPONENTS[rname]['_chem_comp.name']
                })
            chain_info = {
                'id': cid,
                'sequence': seq, 
                'secondary_structure': sst,
                'length': len(ids),
                'residue_ids': ids,
                'uniprot_names': None,
                'uniprot_accession': None,
                'sequence_clusters': {},
                'au_chain_id': meta_data["assembly_chains"][cid],
                'organism': None,
                'cath_homologous_superfamily': None,
                'cath_topology': None,
                'cath_architecture': None,
                'cath_class': None,
                'continuous': True,
                'interacts_with_dna': False,
                "GO_molecular_function": None,
                "GO_biological_process": None,
                "GO_cellular_component": None
            }
            PRO_CHAINS.append(chain_info)
        
        # Add residue map
        RES_MAP = {}
        for R in PRO_RES:
            RES_MAP[R["id"]] = R
        
        # Build Peptides
        PEPTIDES = []
        for chain_info in PRO_CHAINS:
            PEPTIDES += getPeptides(chain_info, model[chain_info["id"]])
        PEPTIDE_MAP = {}
        for peptide in PEPTIDES:
            PEPTIDE_MAP[peptide["id"]] = peptide
        
        # Create SS Elements
        SS_ELEMENTS, HELIX, SHEET = getSSE(PEPTIDES, structure, PRO_CHAINS, RES_MAP)
        
        # Remove Hydrogens from model
        stripH(model)
        
        # Get residue SAP
        SAP = getSAP(model, SAP_CUTOFF_DISTANCE, classifier, REGEXES)
        for s in SAP:
            RES_MAP[s["id"]]["sap_score"] = s["SAP"]
        
        # Add n-mer info. An n-mer is a group of polypeptides which make physical contact
        # with one another in the complex
        opts = freesasa.Parameters(param={'n-slices': 50, 'probe-radius': 1.0})
        modelF = getStructureFromModel(model, classifier=classifier)
        neighbors = nx.Graph()
        
        # Populate neighbors array
        NP = len(PEPTIDES)
        for j in range(NP):
            c1 = PEPTIDES[j]["chain"]
            S1 = getBASA.buildStructure(modelF, PEPTIDES[j]["residue_ids"], in_regex=[[c1]], field_keys=["chain_id"])
            B1 = freesasa.calc(S1, parameters=opts).totalArea()
            #B1 = getSESA(model, PEPTIDES[j]["residue_ids"], classifier, prefix)
            neighbors.add_edge(PEPTIDES[j]["id"], PEPTIDES[j]["id"])
            for k in range(j+1, NP):
                c2 = PEPTIDES[k]["chain"]
                S2 = getBASA.buildStructure(modelF, PEPTIDES[k]["residue_ids"], in_regex=[[c2]], field_keys=["chain_id"])
                B2 = freesasa.calc(S2, parameters=opts).totalArea()
                #B2 = getSESA(model, PEPTIDES[k]["residue_ids"], classifier, prefix)
                S12 = getBASA.buildStructure(modelF, PEPTIDES[j]["residue_ids"]+PEPTIDES[k]["residue_ids"], in_regex=[[c1], [c2]], field_keys=["chain_id"])
                B12 = freesasa.calc(S12, parameters=opts).totalArea()
                #B12 = getSESA(model, PEPTIDES[j]["residue_ids"]+PEPTIDES[k]["residue_ids"], classifier, prefix)
                if((B1+B2-B12)/B12 > 0.02):
                    neighbors.add_edge(PEPTIDES[j]["id"], PEPTIDES[k]["id"])
        entities = []
        ignored = []
        

        nx_entities = list(neighbors.subgraph(c) for c in nx.connected_components(neighbors))
        #for entity in list(nx.connected_component_subgraphs(neighbors)):
        for entity in nx_entities:
            peptides = list(entity.nodes())
            res_ids = []
            for pepid in peptides:
                res_ids += PEPTIDE_MAP[pepid]["residue_ids"]
            if(len(res_ids) < C["MINIMUM_RES_COUNT"]):
                ignored.append(res_ids)
                continue
            entities.append({
                "subunits": entity.number_of_nodes(),
                "segments": peptides,
                "id": getHash(*entity.nodes()),
                "number_of_residues": len(res_ids),
                "residue_ids": res_ids
            })
            
            # Get residue SESA values
            resSESA = getSESA(model, res_ids, classifier, prefix, regexes=REGEXES, total=False)
            for rid in resSESA:
                RES_MAP[rid]["sesa"] = resSESA[rid]
            
            # Get surface residues of each entity
            surface_residues = getSurfaceResidues(res_ids, method='sesa', sesa=resSESA)
            for rid in surface_residues:
                RES_MAP[rid]["surface"] = True
            
            ''' TOO LONG
            # Add CV score to residues
            cvFineScores = computeCV(model, res_ids, surface_residues, CV_FINE_CUTOFF_DISTANCE)
            for cv in cvFineScores:
                RES_MAP[cv["id"]]["cv_fine"] = cv["cv"]
            cvCoarseScores = computeCV(model, res_ids, surface_residues, CV_COARSE_CUTOFF_DISTANCE)
            for cv in cvCoarseScores:
                RES_MAP[cv["id"]]["cv_coarse"] = cv["cv"]
            '''
        # Add data to OUT
        OUT.append({
            'chains': PRO_CHAINS,
            'residues': PRO_RES,
            'entities': entities,
            'secondary_structure_elements': SS_ELEMENTS,
            'segments': PEPTIDES,
            'header_info': {
                'helix_remarks': HELIX,
                'sheet_remarks': SHEET
            }
        })
    
    # Load external data if annotate is true
    if(mmcif_dict is not None):
        addExternalData(OUT, prefix, mmcif_dict, meta_data["assembly_chains"])
    
    # Write output to JSON
    PRO = open("{}-protein.json".format(prefix),'w')
    PRO.write(json.dumps(OUT, indent=None, separators=(',', ':')))
    PRO.close()
    
    # Report secondary structue type
    meta_data["protein_secondary_structure_assignment"] = {
        "method": sse_type
    }
    return OUT, MODEL_DSSP
