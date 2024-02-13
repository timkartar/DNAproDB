import subprocess
import os
import re
import json
import numpy as np
from Bio.PDB import NeighborSearch
from Bio.PDB.vectors import calc_angle
import getBASA
from dnaprodb_utils import CHAIN_RE, RESN_RE, RESI_RE, ATOM_RE
from dnaprodb_utils import residueMoiety
from dnaprodb_utils import nucleotideMoiety
from dnaprodb_utils import log, getHash, getID, getCM, C, roundFloats

VDW_CUTOFF_DISTANCE = C["VDW_CUTOFF_DISTANCE"]
INTERACTION_DISTANCE_CUTOFF = C["INTERACTION_DISTANCE_CUTOFF"]
EFFECTIVE_INTERACTION_ANGLE = C["EFFECTIVE_INTERACTION_ANGLE"]

def getMinDistance(nuc, res):
    mindist = 99999
    minNNDist = 0
    count = 0
    for ra in res:
        if(ra.element == 'H'):
            continue
        amin = 99999
        count += 1
        for na in nuc:
            if(na.element == 'H'):
                continue
            d = na-ra
            mindist = min(d, mindist)
            amin = min(d, amin)
        minNNDist += amin
    
    # residue center of mass
    rCM = getCM(res)
    
    # nucleotide center of mass
    nCM = getCM(nuc)
    
    return float(mindist), float(minNNDist)/count, float(np.linalg.norm(rCM-nCM))

def getInteractingPairs(model, REGEXES):
    atoms = []
    for chain in model.get_list():
        for residue in chain.get_list():
            for atom in residue.get_list(): 
                if(atom.element == 'H'):
                    continue
                atoms.append(atom)
    ns = NeighborSearch(atoms)
    
    nucleotides = set()
    residues = set()
    interactions = {}
    neighbor_pairs = ns.search_all(INTERACTION_DISTANCE_CUTOFF, level='R')
    for n in neighbor_pairs:
        nuc = None
        res = None
        name1 = n[0].get_resname()
        name2 = n[1].get_resname()
        
        # Figure out what this pair is
        if(REGEXES.isDNA(name1)):
            nuc = n[0]
        elif(REGEXES.isProtein(name1)):
            res = n[0]
        if(REGEXES.isDNA(name2)):
            nuc = n[1]
        elif(REGEXES.isProtein(name2)):
            res = n[1]
        
        # Add interaction if DNA-protein
        if(nuc and res):
            nid = getID(residue=nuc)
            rid = getID(residue=res)
            nch, nnum, nins = nid.split('.')
            rch, rnum, rins = rid.split('.')
            nname = nuc.get_resname()
            rname = res.get_resname()
            
            nucleotides.add(nid)
            residues.add(rid)
            iid = nid+'@'+rid
            if(iid in interactions):
                continue
            mindist, meanNNDist, CMDist = getMinDistance(nuc, res)
            interactions[iid] = {
                "res_id": rid,
                "res_name": rname.strip(),
                "res_number": int(rnum),
                "res_chain": rch,
                "res_ins": rins,
                "nuc_id": nid,
                "nuc_name": nname.strip(),
                "nuc_number": int(nnum),
                "nuc_chain": nch,
                "nuc_ins": nins,
                "min_distance": mindist,
                "mean_nn_distance": meanNNDist,
                "cm_distance": CMDist
            }
        else:
            continue
    
    return interactions, list(nucleotides), list(residues)

def getGeometry(pdbid):
    FNULL = open(os.devnull, 'w')
    rc = subprocess.call(['x3dna-snap', '--input={}.pdb'.format(pdbid), '--output={}.snap'.format(pdbid)], stdout=FNULL, stderr=FNULL)
    if(rc == 0 and os.access("{}.snap".format(pdbid), os.R_OK)):
        subprocess.call(['x3dna-snap', '--cleanup'], stdout=FNULL, stderr=FNULL)
        SNPFH = open("{}.snap".format(pdbid),"r").readlines()[2:]
        dre = re.compile('({})\.{}({})(?:\^([A-Z]))?'.format(CHAIN_RE,RESN_RE,RESI_RE))
        hsre = re.compile('List of (\d+) base\/amino-acid pseudo stacks')
        hpre = re.compile('List of (\d+) base/amino-acid pseudo pairs')
        i = 0
        GEOMETRY = []
        while (i < len(SNPFH)):
            hpmatch = hpre.match(SNPFH[i])
            if(hpmatch):
                j = int(hpmatch.group(1));
                for k in range(i+2,i+j+2):
                    fields = SNPFH[k].split()
                    nucm = dre.match(fields[3])
                    resm = dre.match(fields[4])
                    if(nucm and resm):
                        if(nucm.group(3)):
                            nuc_id = getID(nucm.group(1), nucm.group(2), nucm.group(3))
                        else:
                            nuc_id = getID(nucm.group(1), nucm.group(2), ' ')
                        if(resm.group(3)):
                            res_id = getID(resm.group(1), resm.group(2), resm.group(3))
                        else:
                            res_id = getID(resm.group(1), resm.group(2), ' ')
                        GEOMETRY.append({
                            'nuc_id': nuc_id,
                            'res_id': res_id,
                            'geometry': 'pseudo_pair'
                        })
                i = k
                continue
            hsmatch = hsre.match(SNPFH[i])
            if(hsmatch):
                j = int(hsmatch.group(1));
                for k in range(i+2,i+j+2):
                    fields = SNPFH[k].split()
                    nucm = dre.match(fields[3])
                    resm = dre.match(fields[4])
                    if(nucm and resm):
                        if(nucm.group(3)):
                            nuc_id = getID(nucm.group(1), nucm.group(2), nucm.group(3))
                        else:
                            nuc_id = getID(nucm.group(1), nucm.group(2), ' ')
                        if(resm.group(3)):
                            res_id = getID(resm.group(1), resm.group(2), resm.group(3))
                        else:
                            res_id = getID(resm.group(1), resm.group(2), ' ')
                        GEOMETRY.append({
                            'nuc_id': nuc_id,
                            'res_id': res_id,
                            'geometry': 'pseudo_stack'
                        })
                break
            else:
                i += 1
    else:
        log("x3dna-snap failed to run or produce output. Check structure!", pdbid)
    subprocess.call(['x3dna-snap', '--cleanup'], stdout=FNULL, stderr=FNULL)
    if(os.access("{}.snap".format(pdbid), os.R_OK)):
        os.remove("{}.snap".format(pdbid))
    FNULL.close()
    return GEOMETRY

def getEffectiveInteractions(center, neighbors):
    """Computes effective interactions as defined by
    protein dna interface database
    """
    angles = np.zeros((len(neighbors), len(neighbors)))
    cv = center.get_vector()
    vectors = [n.get_vector() for n in neighbors]
    effective_neighbors = []
    for i in range(len(neighbors)):
        for j in range(len(neighbors)):
            angles[i][j] = calc_angle(
                vectors[i],
                cv,
                vectors[j]
            )
        if(angles[i].max() < EFFECTIVE_INTERACTION_ANGLE):
            effective_neighbors.append(neighbors[i])
    
    return effective_neighbors

def getVDW(model, nuc_list, res_list, HBHASH, REGEXES):
    # Get Van der Waal contacts
    atom_list = []
    VDW = []
    for nid in nuc_list:
        chain, nnum, ins = nid.split('.')
        nid = (' ', int(nnum), ins)
        for a in model[chain][nid].get_list():
            if(a.element != 'H'):
                atom_list.append(a)
            
    NAC = len(atom_list) # nucleotide atom count
    for rid in res_list:
        chain, rnum, ins = rid.split('.')
        rid = (' ', int(rnum), ins)
        for a in model[chain][rid].get_list():
            if(a.element != 'H'):
                atom_list.append(a)
    ns = NeighborSearch(atom_list)
    
    for a in atom_list[0:NAC]:
        center = a.get_coord()
        a_id = a.get_full_id()
        a_ch = a_id[2]
        a_resid = a_id[3]
        a_name = a.get_name()
        nuc_id = getID(a_ch, a_resid[1], a_resid[2])
        nucn = model[a_ch][a_resid].get_resname().strip()
        neighbors = ns.search(center, VDW_CUTOFF_DISTANCE)
        res_atoms = []
        for n in neighbors:
            # get a list of residue atoms
            n_id = n.get_full_id()
            n_ch = n_id[2]
            n_resid = n_id[3]
            n_name = n.get_name()
            resn = model[n_ch][n_resid].get_resname().strip()
            if(REGEXES.isProtein(resn)):
                res_atoms.append(n)
        # Include only effective interactions
        effective_neighbors = getEffectiveInteractions(a, res_atoms)
        for n in effective_neighbors:
            n_id = n.get_full_id()
            n_ch = n_id[2]
            n_resid = n_id[3]
            n_name = n.get_name()
            resn = model[n_ch][n_resid].get_resname().strip()
            res_id = getID(n_ch, n_resid[1], n_resid[2])
            hsh = getHash(nuc_id, a_name, res_id, n_name)
            if(hsh not in HBHASH):
                grv = nucleotideMoiety(a_name, nuc_id, REGEXES)
                mty = residueMoiety(n_name, resn, REGEXES)
                VDW.append({
                    'nuc_name': nucn,
                    'nuc_atom': a_name,
                    'nuc_id': nuc_id,
                    'res_name': resn,
                    'res_atom': n_name,
                    'res_id': res_id,
                    'distance': round(a-n,3),
                    'nuc_moiety': grv,
                    'res_moiety': mty
                })
    return VDW

def calculateHBONDS(prefix, DATA_PATH, REGEXES, method="hbplus"):
    # Attempt to run HBPLUS, or fallback to x3dna-snap hbond output.
    FNULL = open(os.devnull, 'w')
    #rc = subprocess.call(['hbadd', '{}.pdb'.format(prefix), os.path.join(DATA_PATH,'components.cif')], stdout=FNULL, stderr=FNULL)
    rc = subprocess.call(['hbplus', '-h', '3.0', '-d', '3.5', '{}.pdb'.format(prefix), '{}.pdb'.format(prefix)], stdout=FNULL, stderr=FNULL)
    FNULL.close()
    HBONDS = {}
    if(rc == 0 and os.access('{}.hb2'.format(prefix), os.R_OK) and method=="hbplus"):
        HB = open('{}.hb2'.format(prefix),'r').readlines()
        for i in range(8,len(HB)):
            d_chain = HB[i][0]
            d_resi = str(int(HB[i][1:5].strip()))
            d_resn = HB[i][6:9].strip()
            d_ins = HB[i][5].replace('-',' ')
            d_atom = HB[i][9:13].strip()
            a_chain = HB[i][14]
            a_resi = str(int(HB[i][15:19].strip()))
            a_ins = HB[i][19].replace('-',' ')
            a_resn = HB[i][20:23].strip()
            a_atom = HB[i][23:27].strip()
            dist = float(HB[i][27:32].strip())
            if(REGEXES.isDNA(d_resn) and REGEXES.isProtein(a_resn)):
                res_id = getID(a_chain, a_resi, a_ins)
                nuc_id = getID(d_chain, d_resi, d_ins)
                grv = nucleotideMoiety(d_atom, nuc_id, REGEXES)
                mty = residueMoiety(a_atom, a_resn, REGEXES)
                key = getHash(nuc_id, d_atom, res_id, a_atom)
                if(key in HBONDS and HBONDS[key]['distance'] <= dist):
                    continue
                HBONDS[key] = {
                    "nuc_atom": d_atom,
                    "nuc_name": d_resn,
                    "res_atom": a_atom,
                    "res_name": a_resn,
                    "distance": dist,
                    "nuc_id": nuc_id,
                    "res_id": res_id,
                    "nuc_moiety": grv,
                    "res_moiety": mty
                }
            elif(REGEXES.isDNA(a_resn) and REGEXES.isProtein(d_resn)):
                res_id = getID(d_chain, d_resi, d_ins)
                nuc_id = getID(a_chain, a_resi, a_ins)
                grv = nucleotideMoiety(a_atom, nuc_id, REGEXES)
                mty = residueMoiety(d_atom, d_resn, REGEXES)
                key = getHash(nuc_id, a_atom, res_id, d_atom)
                if(key in HBONDS and HBONDS[key]['distance'] <= dist):
                    continue
                HBONDS[key] = {
                    "nuc_atom": a_atom,
                    "nuc_name": a_resn,
                    "res_atom": d_atom,
                    "res_name": d_resn,
                    "distance": dist,
                    "nuc_id": nuc_id,
                    "res_id": res_id,
                    "nuc_moiety": grv,
                    "res_moiety": mty
                }
            else:
                continue
        os.remove('{}.hb2'.format(prefix))
    else:
        log('HBPLUS failed to run or produce output. Attempting x3dna-snap for hydrogen bonds.', prefix, Exit=False)
        rc = subprocess.call(['x3dna-snap', '-input={}.pdb'.format(prefix), '-output={}.hbond'.format(prefix), '--get-hbonds'])
        dataRe = re.compile('^({})@({}).({})({})(?:\^([A-Z]))?$'.format(ATOM_RE,CHAIN_RE,RESN_RE,RESI_RE))
        if(rc == 0 and os.access('{}.hbond'.format(prefix), os.R_OK)):
            HB = open('{}.hbond'.format(prefix)).readlines()
            for line in HB[2:]:
                fields = line.split()
                dist = float(fields[4])
                group1 = fields[6]
                group2 = fields[7]
                m1 = dataRe.search(group1)
                m2 = dataRe.search(group2)
                if(m1 and m2):
                    if(REGEXES.isDNA(m1.group(3)) and REGEXES.isProtein(m2.group(3))):
                        if(m1.group(5)):
                            nuc_id = getID(m1.group(2), m1.group(4), m1.group(5))
                        else:
                            nuc_id = getID(m1.group(2), m1.group(4), ' ')
                        if(m2.group(5)):
                            res_id = getID(m2.group(2), m2.group(4), m2.group(5))
                        else:
                            res_id = getID(m2.group(2), m2.group(4), ' ')
                        grv = nucleotideMoiety(m1.group(1), nuc_id, REGEXES)
                        mty = residueMoiety(m2.group(1), m2.group(3), REGEXES)
                        key = getHash(nuc_id, m1.group(1), res_id, m2.group(1))
                        if(key in HBONDS and HBONDS[key]['distance'] <= dist):
                            continue
                        HBONDS[key] = {
                            "nuc_atom": m1.group(1),
                            "nuc_name": m1.group(3),
                            "res_atom": m2.group(1),
                            "res_name": m2.group(3),
                            "distance": dist,
                            "nuc_id": nuc_id,
                            "res_id": res_id,
                            "nuc_moiety": grv,
                            "res_moiety": mty
                        }
                    elif(REGEXES.isDNA(m2.group(3)) and REGEXES.isProteins(m1.group(3))):
                        if(m2.group(5)):
                            nuc_id = getID(m2.group(2), m2.group(4), m2.group(5))
                        else:
                            nuc_id = getID(m2.group(2), m2.group(4), ' ')
                        if(m1.group(5)):
                            res_id = getID(m1.group(2), m1.group(4), m1.group(5))
                        else:
                            res_id = getID(m1.group(2), m1.group(4), ' ')
                        key = getHash(nuc_id, m2.group(1), res_id, m1.group(1))
                        if(key in HBONDS and HBONDS[key]['distance'] <= dist):
                            continue
                        grv = nucleotideMoiety(m2.group(1), nuc_id, REGEXES)
                        mty = residueMoiety(m1.group(1), m1.group(3), REGEXES)
                        HBONDS[key] = {
                            "nuc_atom": m2.group(1),
                            "nuc_name": m2.group(3),
                            "res_atom": m1.group(1),
                            "res_name": m1.group(3),
                            "distance": dist,
                            "nuc_id": nuc_id,
                            "res_id": res_id,
                            "nuc_moiety": grv,
                            "res_moiety": mty
                        }
                    else:
                        log('Residue or Nucleotide not recognized in x3dna-snap hbond output. Check structure file!', prefix, Exit=False)
                        continue
                else:
                    log('Malformed x3dna-snap hbond output. Check if an issue exists!', prefix, Exit=False)
                    continue
            os.remove('{}.hbond'.format(prefix))
        else:
            log('HBPLUS and x3dna-snap failed to run.', prefix)
    return list(HBONDS.values())

def splitEnsemble(prefix, N, REGEXES):
    """ Docstring """
    i = 0
    if(N == 1):
        cName = "{}".format(prefix)
        dName = "{}-DNA".format(prefix)
        pName = "{}-protein".format(prefix)
        hName = "{}-noH".format(prefix)
        
        yield cName, dName, pName, hName, i
    else:
        START_RE = re.compile('^MODEL')
        STOP_RE = re.compile('^ENDMDL')
        FH = open("{}.pdb".format(prefix))
        cName = "{}_{}".format(prefix, i)
        dName = "{}_{}-DNA".format(prefix, i)
        pName = "{}_{}-protein".format(prefix, i)
        hName = "{}_{}-noH".format(prefix, i)
        
        FOUT = open("{}.pdb".format(cName), "w")
        DOUT = open("{}.pdb".format(dName), "w")
        POUT = open("{}.pdb".format(pName), "w")
        HOUT = open("{}.pdb".format(hName), "w")
        for line in FH:
            if(START_RE.search(line)):
                continue
            elif(STOP_RE.search(line)):
                FOUT.close()
                DOUT.close()
                POUT.close()
                HOUT.close()
                yield cName, dName, pName, hName, i
                os.remove("{}.pdb".format(cName))
                os.remove("{}.pdb".format(dName))
                os.remove("{}.pdb".format(pName))
                os.remove("{}.pdb".format(hName))
                
                i += 1
                if(i == N):
                    break
                cName = "{}_{}".format(prefix, i)
                dName = "{}_{}-DNA".format(prefix, i)
                pName = "{}_{}-protein".format(prefix, i)
                hName = "{}_{}-noH".format(prefix, i)
                FOUT = open("{}.pdb".format(cName), "w")
                DOUT = open("{}.pdb".format(dName), "w")
                POUT = open("{}.pdb".format(pName), "w")
                HOUT = open("{}.pdb".format(hName), "w")
            else:
                FOUT.write(line)
                res = line[17:20].strip()
                elm = line[76:78].strip()
                if(REGEXES.isDNA(res)):
                    DOUT.write(line)
                elif(REGEXES.isProtein(res)):
                    POUT.write(line)
                if(elm != 'H'):
                    HOUT.write(line)
        FH.close()
    # Clean up left over junk
    remove = [
        "bp_order.dat",
        "auxiliary.par",
        "bp_helical.par",
        "bp_step.par",
        "bestpairs.pdb",
        "hbadd.bonds",
        "stacking.pdb",
        "ref_frames.dat",
        "poc_haxis.r3d",
        "hstacking.pdb",
        "hel_regions.pdb",
        "hbdebug.dat",
        "cf_7methods.par",
        "hbadd.map",
        "hbadd.sum",
        "hbplus.rc"
    ]
    for r in remove:
        if(os.access(r, os.R_OK)):
            os.remove(r)

def process(prefix, N, COMPONENTS, assembly, DSSP, DATA_PATH, REGEXES, NUCLEOTIDES, IDS):
    # Generate ID array
    OUT = []
    
    # Loop over every model in ensemble
    for c,d,p,h,i in splitEnsemble(prefix, N, REGEXES):        
        REGEXES.setModel(i)
        interactions = {
            "basa": None,
            "hbond": None,
            "vdw": None,
            "geometry": None,
            "nucleotide-residue_interactions": None
        }
        
        # Get pair list
        int_pairs, nuc_list, res_list = getInteractingPairs(assembly[i], REGEXES)
        if(len(int_pairs) == 0):
            log("No nucleotide-residue pairs meet the interaction cut-off threshold.", prefix)
        interactions["nucleotide-residue_interactions"] = list(int_pairs.values())
        interface_ids = {
            "res_ids": res_list,
            "nuc_ids": nuc_list,
            "complex_ids": res_list + nuc_list
        }
        
        ## Determine a list of DNA and protein entity pairs for consideration
        #dna_entity_lookup = {}
        #pro_entity_lookup = {}
        #pairs = set()
        #for nr in interactions["nucleotide-residue_interactions"]:
            #pairs.add((dna_entity_lookup[nr["nuc_id"]], pro_entity_lookup[nr["pro_id"]]))
        
        # Perform BASA calculations
        interactions['basa'] = getBASA.basa(assembly[i], COMPONENTS, REGEXES, DATA_PATH, NUCLEOTIDES[i], IDS, interface_ids, dssp=DSSP[i])
        
        # Get Hydrogen Bonds and VdW contacts
        HBONDS = calculateHBONDS(h, DATA_PATH, REGEXES)
        VDW = getVDW(assembly[i], nuc_list, res_list, HBONDS, REGEXES)
        interactions['hbond'] = HBONDS
        interactions['vdw'] = VDW
        
        # Get Residue-Nucleotide interaction geometry
        GEO = getGeometry(c)
        interactions['geometry'] = GEO
        
        OUT.append(interactions)
    
    # Write data to file
    for o in OUT:
        roundFloats(o)
    IOUT = open("{}-interactions.json".format(prefix),"w")
    IOUT.write(json.dumps(OUT,indent=None,separators=(',', ':'),sort_keys=True))
    IOUT.close()
    return OUT
