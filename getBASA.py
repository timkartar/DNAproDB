import freesasa
import re
import os
import subprocess
import json
import math
from Bio.PDB import NeighborSearch
from dnaprodb_utils import CHAIN_RE, RESN_RE, RESI_RE, ATOM_RE, ID_RE, INS_RE
from dnaprodb_utils import roundFloats
from dnaprodb_utils import residueMoiety
from dnaprodb_utils import nucleotideMoiety
from dnaprodb_utils import getCM
from dnaprodb_utils import log, getID, C
from dnaprodb_utils import getStructureFromModel

nucGrvLabel = C["NUC_MTY_LABEL_DS"]
nucMtyLabel = C["NUC_MTY_LABEL_SS"]
resMtyLabel = C["RES_MTY_LABEL"]
resSST = C["RES_SST"]
regexFieldKeyMap = {
    "atom_name": re.compile(ATOM_RE),
    "res_name": re.compile(RESN_RE), 
    "res_num": re.compile(RESI_RE), 
    "chain_id": re.compile(CHAIN_RE),
    "res_id": re.compile(ID_RE),
    "ins_code": re.compile(INS_RE)
}

class MatchField(object):
    def __init__(self, components=None, nucleotides=None, field=None, field_values=None):
        self.flag = None # set to decide what to do with search
        if(components is not None):
            if(field is None):
                raise ValueError("Argument 'field' must be set when \
                specifying component look-up!")
            self.components = components
            self.flag = "component_lookup"
            self.field = field
            self.field_value = field_value
        elif(nucleotides is not None):
            if(field is None):
                raise ValueError("Argument 'field' must be set when \
                specifying nucleotide lookup!")
            self.nucleotides = nucleotides
            self.flag = "nucleotide_lookup"
            
    def search(self, arg):
        if(self.flag == "component_lookup"):
            return self.components[arg][self.field] in self.field_values
        elif(self.flag == "nucleotide_lookup"):
            return self.nucleotides[arg][self.field] in self.field_values

class GeneralClassifier(freesasa.Classifier):
    def initialize(self, DATA_PATH, COMPONENTS, fileName='vdw-radii.json'):
        self.components = COMPONENTS
        with open(os.path.join(DATA_PATH, fileName)) as FILE:
            self.radii = json.load(FILE)
    
    def radius(self, residueName, atomName):
        rName = residueName.strip()
        aName = atomName.strip()
        if(rName in self.radii):
            # Standard Residue
            if(aName in self.radii[rName]):
                return self.radii[rName][aName]
            else:
                return self.radii['element'][self.getElement(rName, atomName)]
        elif(rName in self.components):
            # Non-standard known residue
            parent = self.components[rName]['_chem_comp.mon_nstd_parent_comp_id']
            if(parent in self.radii and aName in self.radii[parent]):
                return self.radii[parent][aName]
            else:
                return self.radii['element'][self.getElement(rName, atomName)]
        else:
            # Unknown residue - make best guess for atom element
            print(("Unknown residue: {}".format(rName)))
            return self.radii['element'][self.guessElement(atomName)]
    
    def classify(self, residueName, atomName):
        return "atom"
    
    def getElement(self, residueName, atomName):
        aName = atomName.strip()
        if(residueName in self.components):
            try:
                index = self.components[residueName]['_chem_comp_atom.atom_id'].index(aName)
                return self.components[residueName]['_chem_comp_atom.type_symbol'][index]
            except:
                return self.guessElement(atomName)
        else:
            return self.guessElement(atomName)
    
    def guessElement(self, atomName):
        """Tries to guess element from atom name if not recognised."""
        print(("Got :{}".format(atomName)))
        name = atomName.strip()
        if name.capitalize() not in self.radii["element"]:
            # Inorganic elements have their name shifted left by one position
            #  (is a convention in PDB, but not part of the standard).
            # isdigit() check on last two characters to avoid mis-assignment of
            # hydrogens atoms (GLN HE21 for example)
            if atomName[0].isalpha() and not (atomName[2:].isdigit() or atomName[2:] == "''"):
                putative_element = name
            else:
                # Hs may have digit in first position
                if name[0].isdigit():
                    putative_element = name[1]
                else:
                    putative_element = name[0]
        
            if putative_element.capitalize() in self.radii["element"]:
                element = putative_element
            else:
                element = ""
            
            print(("Guessed(e): {}".format(element)))
            return element
        else:
            print(("Guessed(n): {}".format(name)))
            return name

def sumSASA(structure, result, REGEXES, IDs):
    N = structure.nAtoms()
    SUM = {}
    for i in range(N):
        resn = structure.residueName(i).strip()
        resi = structure.residueNumber(i).strip()
        chain = structure.chainLabel(i)
        aname = structure.atomName(i).strip()
        if(resi[-1].isdigit()):
            ins = " "
        else:
            ins = resi[-1]
            resi = resi[:-1]
        rid = getID(chain, resi, ins)
        if(rid not in IDs):
            continue
        
        # Get regexes
        if(REGEXES.isProtein(resn)):
            # This is a peptide-linking residue
            classes = REGEXES[resn]
            nc = len(classes)
        elif(REGEXES.isDNA(resn)):
            # This is a DNA-linking nucleotide
            classes = REGEXES[rid]
            nc = len(classes)
        else:
            raise ValueError("Neither residue name '{}' or residue id '{}' found!".format(resn, rid))
        sasa = result.atomArea(i)
        if(rid not in SUM):
            SUM[rid] = {
                'sasa': [0.0]*(nc+1),
                'resn': resn,
                'chain': chain,
                'resi': resi,
                'ins': ins,
                'id': rid
            }
        SUM[rid]['sasa'][0] += sasa
        for j in range(nc):
            if(classes[j].search(aname)):
                SUM[rid]['sasa'][j+1] += sasa
                break
    return SUM

def makeRegexField(field_key, aname, resn, resi, chain, ins):
    if(field_key == "atom_name"):
        return aname 
    elif(field_key == "res_name"):
        return resn
    elif(field_key == "res_num"):
        return resi
    elif(field_key == "chain_id"):
        return chain
    elif(field_key == "res_id"):
        return "{}.{}.{}".format(chain, resi, ins)
    elif(field_key == "ins_code"):
        return ins

def buildStructure(template, IDs,
    in_regex=None,
    ex_regex=None, 
    field_keys=["atom_name", "res_name", "res_num", "chain_id"]
    ):
    if(in_regex is None):
        in_regex = [[]]
        for fk in field_keys:
            in_regex[0].append(regexFieldKeyMap[fk])
    
    if(ex_regex is None):
        ex_regex = [['(?!.*)']*len(field_keys)]
    
    # Compile Expressions
    for regex in in_regex:
        # check that regex list matches field_keys
        if(len(regex) != len(field_keys)):
            raise ValueError("Length of regex lists and field_keys must match!")
        for i in range(len(regex)):
            if(isinstance(regex[i], re.Pattern) or isinstance(regex[i], MatchField)):
                continue
            regex[i] = re.compile(regex[i])
    
    for regex in ex_regex:
        # check that regex list matches field_keys
        if(len(regex) != len(field_keys)):
            raise ValueError("Length of regex lists and field_keys must match!")
        for i in range(len(regex)):
            if(isinstance(regex[i], re.Pattern) or isinstance(regex[i], MatchField)):
                continue
            regex[i] = re.compile(regex[i])
    
    N = template.nAtoms()
    structure = freesasa.Structure()
    for i in range(N):
        IMATCH = False
        EMATCH = False
        resn = template.residueName(i)
        resi = template.residueNumber(i)
        chain = template.chainLabel(i)
        aname = template.atomName(i)
        resn_s = resn.strip()
        resi_s = resi.strip()
        aname_s = aname.strip()
        if(resi_s[-1].isdigit()):
            ins = " "
        else:
            ins = resi_s[-1]
            resi_s = resi_s[:-1]
        rid = getID(chain, resi_s, ins)
        if(rid not in IDs):
            continue
        
        for regex in in_regex:
            matches = []
            for j in range(len(regex)):
                # every regex will be the same length as field_keys
                matches.append(regex[j].search(makeRegexField(field_keys[j], aname_s, resn_s, resi_s, chain, ins)))
            if all(matches):
                IMATCH = True
                break
        
        for regex in ex_regex:
            matches = []
            for j in range(len(regex)):
                # every regex will be the same length as field_keys
                matches.append(regex[j].search(makeRegexField(field_keys[j], aname_s, resn_s, resi_s, chain, ins)))
            if all(matches):
                EMATCH = True
                break
        
        if(IMATCH and not EMATCH):
            coord = template.coord(i)
            structure.addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
    
    return structure

def getComplexBASA(model, classifier, REGEXES, NUCLEOTIDES, IDS, INT_IDS, dssp=None, detailed=False):
    # Initialize Structures
    opts = freesasa.Parameters(param={'n-slices': 50})
    try:
        com = getStructureFromModel(model, classifier)
    except AssertionError:
        log("Invalid atom/residue combination. All atom/residue names should follow standard PDB conventions.", model.get_parent().get_id())
    
    
    pro = freesasa.Structure()
    dna = freesasa.Structure()
    if(detailed):
        com_helix = freesasa.Structure() # complex with only helical DNA
        com_other = freesasa.Structure() # complex with only non-helical DNA
        # Helical DNA-protein complex
        pro_grc = [ # [ wg, sg, sr, pp ]
            freesasa.Structure(),
            freesasa.Structure(),
            freesasa.Structure(),
            freesasa.Structure()
        ]
        
        # Non-helical DNA-protein complex 
        pro_ssc = [ # [ bs, sr, pp ]
            freesasa.Structure(),
            freesasa.Structure(),
            freesasa.Structure()
        ]
        
        # DNA-SSE complex
        dna_ssc = [ # [ H, S, L ]
            freesasa.Structure(),
            freesasa.Structure(),
            freesasa.Structure()
        ]
    
    # Generate ids of all nearby residues and nucleotides in the interface
    # (this increases perfomance)
    COMPLEX_IDS = set()
    atoms = []
    for chain in model.get_list():
        for residue in chain.get_list():
            for atom in residue.get_list(): 
                if(atom.element == 'H'):
                    continue
                atoms.append(atom)
    ns = NeighborSearch(atoms)
    
    for resID in INT_IDS["complex_ids"]:
        ch, resi, ins = resID.split('.')
        residue = model[ch][(' ', int(resi), ins)]
        cm = getCM(residue)
        neighbors = ns.search(cm, 20.0, level='R')
        for n in neighbors:
            nid = getID(residue=n)
            COMPLEX_IDS.add(nid)
    
    # Build Structures
    N = com.nAtoms()
    for i in range(N):
        resn = com.residueName(i)
        if(REGEXES.isProtein(resn)):
            coord = com.coord(i)
            resi = com.residueNumber(i)
            chain = com.chainLabel(i)
            aname = com.atomName(i)
            
            # Add atoms to relevant structures
            pro.addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
            if(detailed and rid in COMPLEX_IDS):
                resi_s = resi.strip()
                if(resi_s[-1].isdigit()):
                    ins = " "
                else:
                    ins = resi_s[-1]
                    resi_s = resi_s[:-1]
                rid = getID(chain, resi_s, ins)
                if(rid in dssp):
                    ss = dssp[rid]
                else:
                    ss = 'L'
                com_helix.addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                com_other.addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                for j in range(len(pro_grc)):
                    pro_grc[j].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                for j in range(len(pro_ssc)):
                    pro_ssc[j].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                if(ss == 'H'):
                    # exclude helix atoms
                    dna_ssc[1].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                    dna_ssc[2].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                elif(ss == 'S'):
                    # exclude strand atoms
                    dna_ssc[2].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                    dna_ssc[0].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                elif(ss == 'L'):
                    # exclude loop atoms
                    dna_ssc[0].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                    dna_ssc[1].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
        elif(REGEXES.isDNA(resn)):
            coord = com.coord(i)
            resi = com.residueNumber(i)
            chain = com.chainLabel(i)
            aname = com.atomName(i)
            
            # Add atoms to relevant structures
            dna.addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
            if(detailed and nid in COMPLEX_IDS):
                aname_s = aname.strip()
                resi_s = resi.strip()
                if(resi_s[-1].isdigit()):
                    ins = " "
                else:
                    ins = resi_s[-1]
                    resi_s = resi_s[:-1]
                nid = getID(chain, resi_s, ins)
                for j in range(3):
                    dna_ssc[j].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                if(nid in NUCLEOTIDES and NUCLEOTIDES[nid]["secondary_structure"] == "helical"):
                    com_helix.addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                    # Helical DNA (Major/Minor groove defined)
                    if(REGEXES[nid][0].search(aname_s)):
                        # exclude Major groove atoms
                        pro_grc[1].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                        pro_grc[2].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                        pro_grc[3].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                    elif(REGEXES[nid][1].search(aname_s)):
                        # exclude Minor groove atoms
                        pro_grc[0].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                        pro_grc[2].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                        pro_grc[3].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                    elif(REGEXES[nid][2].search(aname_s)):
                        # exclude Sugar atoms
                        pro_grc[0].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                        pro_grc[1].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                        pro_grc[3].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                    elif(REGEXES[nid][3].search(aname_s)):
                        # exclude Phosphate atoms
                        pro_grc[0].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                        pro_grc[1].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                        pro_grc[2].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                else:
                    com_other.addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                    # Non-helical DNA
                    if(REGEXES[nid][0].search(aname_s)):
                        # exclude Base atoms
                        pro_ssc[1].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                        pro_ssc[2].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                    elif(REGEXES[nid][1].search(aname_s)):
                        # exclude Sugar atoms
                        pro_ssc[0].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                        pro_ssc[2].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                    elif(REGEXES[nid][2].search(aname_s)):
                        # exclude Phosphate atoms
                        pro_ssc[0].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
                        pro_ssc[1].addAtom(aname, resn, resi, chain, coord[0], coord[1], coord[2])
    
    # Set up classifiers for each structure
    com.setRadiiWithClassifier(classifier)
    pro.setRadiiWithClassifier(classifier)
    dna.setRadiiWithClassifier(classifier)
    if(detailed):
        com_helix.setRadiiWithClassifier(classifier)
        com_other.setRadiiWithClassifier(classifier)
        for i in range(len(pro_grc)):
            pro_grc[i].setRadiiWithClassifier(classifier)
        for i in range(len(pro_ssc)):
            pro_ssc[i].setRadiiWithClassifier(classifier)
        for i in range(len(dna_ssc)):
            dna_ssc[i].setRadiiWithClassifier(classifier)
    
    # Run BASA calculations
    com_sasa = freesasa.calc(com, parameters=opts)
    pro_sasa = freesasa.calc(pro, parameters=opts)
    dna_sasa = freesasa.calc(dna, parameters=opts)
    if(detailed):
        com_helix_sasa = freesasa.calc(com_helix, parameters=opts)
        com_other_sasa = freesasa.calc(com_other, parameters=opts)
        dna_ssc_sasa = []
        pro_grc_sasa = []
        pro_ssc_sasa = []
        for i in range(len(pro_grc)):
            pro_grc_sasa.append(freesasa.calc(pro_grc[i], parameters=opts))
        for i in range(len(pro_ssc)):
            pro_ssc_sasa.append(freesasa.calc(pro_ssc[i], parameters=opts))
        for i in range(len(dna_ssc)):
            dna_ssc_sasa.append(freesasa.calc(dna_ssc[i], parameters=opts))
    
    IN_LIST = []
    # Get Residue BASA values
    RES = {}
    RES_COMPLEX_SASA = sumSASA(com, com_sasa, REGEXES, IDS['protein'])
    RES_PROTEIN_SASA = sumSASA(pro, pro_sasa, REGEXES, IDS['protein'])
    if(detailed):
        RES_COM_HELIX_SASA = sumSASA(com_helix, com_helix_sasa, REGEXES, IDS['protein'])
        RES_COM_OTHER_SASA = sumSASA(com_other, com_other_sasa, REGEXES, IDS['protein'])
        RES_GRC_SASA = []
        RES_SSC_SASA = []
        for i in range(len(pro_grc)):
            RES_GRC_SASA.append(sumSASA(pro_grc[i], pro_grc_sasa[i], REGEXES, IDS['protein']))
        for i in range(len(pro_ssc)):
            RES_SSC_SASA.append(sumSASA(pro_ssc[i], pro_ssc_sasa[i], REGEXES, IDS['protein']))
    
    for key in RES_COMPLEX_SASA:
        basa = RES_PROTEIN_SASA[key]['sasa'][0]-RES_COMPLEX_SASA[key]['sasa'][0]
        # Compute the various BASA, SASA and FASA components
        RES[key] = {
            'resn': RES_COMPLEX_SASA[key]['resn'],
            'resi': RES_COMPLEX_SASA[key]['resi'],
            'id': RES_COMPLEX_SASA[key]['id'],
            'ins': RES_COMPLEX_SASA[key]['ins'],
            'chain': RES_COMPLEX_SASA[key]['chain'],
            'fasa': {
                'total': RES_PROTEIN_SASA[key]['sasa'][0],
                'mc': RES_PROTEIN_SASA[key]['sasa'][1],
                'sc': RES_PROTEIN_SASA[key]['sasa'][2]
            },
            'sasa':{
                'total': RES_COMPLEX_SASA[key]['sasa'][0],
                'mc': RES_COMPLEX_SASA[key]['sasa'][1],
                'sc':RES_COMPLEX_SASA[key]['sasa'][2]
            },
            'basa': { 
                'total': basa,
                'mc': RES_PROTEIN_SASA[key]['sasa'][1]-RES_COMPLEX_SASA[key]['sasa'][1],
                'sc': RES_PROTEIN_SASA[key]['sasa'][2]-RES_COMPLEX_SASA[key]['sasa'][2],
            }
        }
        
        # Set up keys and values
        if(detailed):
            RES[key]['basa']['dna_moieties'] = {}
            for i in range(len(resMtyLabel)):
                RES[key]['basa']['dna_moieties'][resMtyLabel[i]] = {}
                for j in range(len(nucGrvLabel)):
                    RES[key]['basa']['dna_moieties'][resMtyLabel[i]][nucGrvLabel[j]] = 0.0
                for j in range(len(nucMtyLabel)):
                    RES[key]['basa']['dna_moieties'][resMtyLabel[i]][nucMtyLabel[j]] = 0.0
            
            # Helical DNA segements
            for i in range(len(resMtyLabel)):
                for j in range(len(nucGrvLabel)):
                    if(key in INT_IDS["res_ids"]):
                        RES[key]['basa']['dna_moieties'][resMtyLabel[i]][nucGrvLabel[j]] += RES_GRC_SASA[j][key]['sasa'][i+1]-RES_COM_HELIX_SASA[key]['sasa'][i+1]
            
            # Non-helical DNA segements
            for i in range(len(resMtyLabel)):
                for j in range(len(nucMtyLabel)):
                    if(key in INT_IDS["res_ids"]):
                        RES[key]['basa']['dna_moieties'][resMtyLabel[i]][nucMtyLabel[j]] += RES_SSC_SASA[j][key]['sasa'][i+1]-RES_COM_OTHER_SASA[key]['sasa'][i+1]
    
    # Get Nucleotide BASA values
    NUC = {}
    NUC_COMPLEX_SASA = sumSASA(com, com_sasa, REGEXES, IDS['dna'])
    NUC_DNA_SASA = sumSASA(dna, dna_sasa, REGEXES, IDS['dna'])
    if(detailed):
        NUC_SSC_SASA = []
        for i in range(3):
            NUC_SSC_SASA.append(sumSASA(dna_ssc[i], dna_ssc_sasa[i], REGEXES, IDS['dna']))
    
    for key in NUC_COMPLEX_SASA:
        basa = NUC_DNA_SASA[key]['sasa'][0]-NUC_COMPLEX_SASA[key]['sasa'][0]
        NUC[key] = {
            'nucn': NUC_COMPLEX_SASA[key]['resn'],
            'nuci': NUC_COMPLEX_SASA[key]['resi'],
            'chain': NUC_COMPLEX_SASA[key]['chain'],
            'id': NUC_COMPLEX_SASA[key]['id'],
            'ins': NUC_COMPLEX_SASA[key]['ins'],
            'fasa': {
                'total': NUC_DNA_SASA[key]['sasa'][0],
            },
            'sasa': {
                'total': NUC_COMPLEX_SASA[key]['sasa'][0],
            },
            'basa': {
                'total': basa
            }
        }
        
        if(NUCLEOTIDES[key]["secondary_structure"] == "helical"):
            labels = nucGrvLabel
        else:
            labels = nucMtyLabel
        
        for i in range(1, len(labels)+1):
            NUC[key]['fasa'][labels[i-1]] = NUC_DNA_SASA[key]['sasa'][i]
            NUC[key]['sasa'][labels[i-1]] = NUC_COMPLEX_SASA[key]['sasa'][i]
            NUC[key]['basa'][labels[i-1]] = NUC_DNA_SASA[key]['sasa'][i] - NUC_COMPLEX_SASA[key]['sasa'][i]
        
        if(detailed):
            NUC[key]['basa']['secondary_structure'] = {}
            for i in range(1, len(labels)+1):
                NUC[key]['basa']['secondary_structure'][labels[i-1]] = {}
                for j in range(len(resSST)):
                    if(key in INT_IDS["nuc_ids"]):
                        NUC[key]['basa']['secondary_structure'][labels[i-1]][resSST[j]] = NUC_SSC_SASA[j][key]['sasa'][i]-NUC_COMPLEX_SASA[key]['sasa'][i]
                    else:
                        NUC[key]['basa']['secondary_structure'][labels[i-1]][resSST[j]] = 0.0
    
    # Calculate Interaction BASA
    INT = {}
    # Calculate the components B^(Nj)_g(Ri) and B^(Nj)_g,m(Ri), the BASA components of nucleotide Nj due to residue Ri
    for rkey in RES:
        if(RES[rkey]['basa']['total'] > 0 and rkey in INT_IDS["res_ids"]):
            # Calculate B^(Nj)_g(Ri)
            S = buildStructure(com, COMPLEX_IDS, ex_regex=[[rkey]], field_keys=["res_id"])
            S.setRadiiWithClassifier(classifier)
            SASA = sumSASA(S, freesasa.calc(S, parameters=opts), REGEXES, IDS['dna'])
            for nkey in SASA:
                if(nkey in NUC and nkey in INT_IDS["nuc_ids"] and NUC[nkey]['basa']['total'] > 0):
                    if(NUCLEOTIDES[nkey]["secondary_structure"] == "helical"):
                        nucLabels = nucGrvLabel
                    else:
                        nucLabels = nucMtyLabel
                    INT[rkey+nkey] = {
                        'res_id': rkey,
                        'res_name': RES[rkey]['resn'],
                        'res_number': RES[rkey]['resi'],
                        'res_chain': RES[rkey]['chain'],
                        'nuc_id': nkey,
                        'nuc_name': NUC[nkey]['nucn'],
                        'nuc_number': NUC[nkey]['nuci'],
                        'nuc_chain': NUC[nkey]['chain'],
                        'nuc_basa': {
                            'total': SASA[nkey]['sasa'][0] - NUC[nkey]['sasa']['total'],
                        },
                        'res_basa': {
                            'total': None,
                            'mc': None,
                            'sc': None
                        },
                        'basa': {
                            'total': SASA[nkey]['sasa'][0] - NUC[nkey]['sasa']['total'],
                        }
                    }
                    
                    # Add nucleotide BASA
                    for i in range(1, len(nucLabels)+1):
                        INT[rkey+nkey]['nuc_basa'][nucLabels[i-1]] = SASA[nkey]['sasa'][i] - NUC[nkey]['sasa'][nucLabels[i-1]]
                    
                    # Add joint BASA
                    for i in range(len(nucLabels)):
                        INT[rkey+nkey]['basa'][nucLabels[i]] = {}
                        for j in range(len(resMtyLabel)):
                            INT[rkey+nkey]['basa'][nucLabels[i]][resMtyLabel[j]] = 0.0
            
            # Calculate B^(Nj)_g,m(Ri)
            for i in range(len(resMtyLabel)):
                m = resMtyLabel[i]
                S = buildStructure(com, COMPLEX_IDS, ex_regex=[[REGEXES[RES[rkey]['resn']][i], rkey]], field_keys=["atom_name", "res_id"])
                S.setRadiiWithClassifier(classifier)
                SASA = sumSASA(S, freesasa.calc(S, parameters=opts), REGEXES, IDS['dna'])
                for nkey in SASA:
                    if(nkey in NUC and NUC[nkey]['basa']['total'] > 0 and nkey in INT_IDS["nuc_ids"]):
                        if(NUCLEOTIDES[nkey]["secondary_structure"] == "helical"):
                            nucLabels = nucGrvLabel
                        else:
                            nucLabels = nucMtyLabel
                        for j in range(len(nucLabels)):
                            g = nucLabels[j]
                            INT[rkey+nkey]['basa'][g][m] += SASA[nkey]['sasa'][j+1] - NUC[nkey]['sasa'][g]
    
    # Calculate the components B^(Ri)_m(Nj) and B^(Ri)_m,g(Nj), the BASA components of residue Ri due to nucleotide Nj
    for nkey in NUC:
        if(nkey in INT_IDS["nuc_ids"] and NUC[nkey]['basa']['total'] > 0):
            if(NUCLEOTIDES[nkey]["secondary_structure"] == "helical"):
                nucLabels = nucGrvLabel
            else:
                nucLabels = nucMtyLabel
            S = buildStructure(com, COMPLEX_IDS, ex_regex=[[nkey]], field_keys=["res_id"])
            S.setRadiiWithClassifier(classifier)
            SASA = sumSASA(S, freesasa.calc(S, parameters=opts), REGEXES, IDS['protein'])
            for rkey in SASA:
                if(rkey in RES and rkey in INT_IDS["res_ids"] and RES[rkey]['basa']['total'] > 0):
                    INT[rkey+nkey]['res_basa']['total'] = SASA[rkey]['sasa'][0] - RES[rkey]['sasa']['total']
                    INT[rkey+nkey]['res_basa']['mc'] = SASA[rkey]['sasa'][1] - RES[rkey]['sasa']['mc']
                    INT[rkey+nkey]['res_basa']['sc'] = SASA[rkey]['sasa'][2] - RES[rkey]['sasa']['sc']
                    INT[rkey+nkey]['basa']['total'] += SASA[rkey]['sasa'][0] - RES[rkey]['sasa']['total']
            
            # Calculate B^(Ri)_m,g(Nj)
            for i in range(len(nucLabels)):
                g = nucLabels[i]
                S = buildStructure(com, COMPLEX_IDS, ex_regex=[[REGEXES[nkey][i], nkey]], field_keys=["atom_name", "res_id"])
                S.setRadiiWithClassifier(classifier)
                SASA = sumSASA(S, freesasa.calc(S, parameters=opts), REGEXES, IDS['protein'])
                for rkey in SASA:
                    if(rkey in RES and RES[rkey]['basa']['total'] > 0 and rkey in INT_IDS["res_ids"]):
                        for j in range(len(resMtyLabel)):
                            m = resMtyLabel[j]
                            INT[rkey+nkey]['basa'][g][m] += SASA[rkey]['sasa'][j+1] - RES[rkey]['sasa'][m]
    
    # Remove indirect interactions (i.e. where one or both have no non-overlapping BASA components
    deleteList = []
    for key in INT:
        if(INT[key]['res_basa']['total'] <= 0 or INT[key]['nuc_basa']['total'] <= 0):
            deleteList.append(key)
    for key in deleteList:
        del INT[key]
    
    BASA = {
        'residues': list(RES.values()),
        'nucleotides': list(NUC.values()),
        'interactions': list(INT.values())
    }
    roundFloats(BASA,3)
    
    return BASA

def basa(model, COMPONENTS, REGEXES, DATA_PATH, NUCLEOTIDES, IDS, INT_IDS, dssp=None, quiet=True, detailed=False):
    if(quiet):
        freesasa.setVerbosity(freesasa.nowarnings)
    classifier = GeneralClassifier()
    classifier.initialize(DATA_PATH, COMPONENTS)
    
    BASA = getComplexBASA(model, classifier, REGEXES, NUCLEOTIDES, IDS, INT_IDS, dssp=dssp, detailed=detailed)
    return BASA

#### DEBUGGING CODE ###
#from Bio.PDB import PDBIO
#from Bio.PDB.Structure import Structure
#from Bio.PDB.Model import Model
#from Bio.PDB.Chain import Chain
#from Bio.PDB.Residue import Residue
#from Bio.PDB.Atom import Atom
#from Bio.Data import IUPACData
#def freeSASA2PDB(structure):
    #OUT = Structure('out')
    #OUT.add(Model(0))
    #N = structure.nAtoms()
    #for i in range(N):
        #resn = structure.residueName(i).strip()
        #resi = structure.residueNumber(i).strip()
        #chain = structure.chainLabel(i)
        #aname = structure.atomName(i).strip()
        #coord = structure.coord(i)
        #element = assignElement(aname)
        #if(len(aname) < 4):
            #am = re.search("{}([A-Z0-9'\*]+)?".format(element), aname)
            #if(aname[-1] == "'"):
                #fname = " {}".format(aname)
            #elif(am.group(1)):
                #fname = "{:>2s}{:<s}".format(element, am.group(1))
            #else:
                #fname = "{:>2s}  ".format(element)
        #else:
            #aname = fname
        
        ## check if chain there
        #if(not (chain in OUT[0])):
            #OUT[0].add(Chain(chain))
        
        ## check if residue there
        #rid = (' ', int(resi), ' ')
        #if(not (rid in OUT[0][chain])):
            #OUT[0][chain].add(Residue(rid, resn, resi))
        
        ## add atom
        #OUT[0][chain][rid].add(Atom(aname, coord, 1.0, 1.0, ' ', fname, i, element=element))
    
    ## write OUT to file
    #io = PDBIO()
    ##io.set_structure(OUT)
    ##io.save("freesasa.pdb")

#def assignElement(fullname):
    #"""Tries to guess element from atom name if not recognised."""
    #name = fullname.strip()
    #if name.capitalize() not in IUPACData.atom_weights:
        ## Inorganic elements have their name shifted left by one position
        ##  (is a convention in PDB, but not part of the standard).
        ## isdigit() check on last two characters to avoid mis-assignment of
        ## hydrogens atoms (GLN HE21 for example)
    
        #if fullname[0].isalpha() and not (fullname[2:].isdigit() or fullname[2:] == "''"):
            #putative_element = name
        #else:
            ## Hs may have digit in first position
            #if name[0].isdigit():
                #putative_element = name[1]
            #else:
                #putative_element = name[0]
    
        #if putative_element.capitalize() in IUPACData.atom_weights:
            #element = putative_element
        #else:
            #element = ""
        
        #return element
    #else:
        #return name[0]
#### DEBUGGING CODE ####

# Removed Code - calculate direct RES-NUC basa
#for rkey in RES:
    #S = buildStructure(com, 
        #in_regex=[
            #[ATM, RES[rkey]['resn'], RES[rkey]['resi'] ,RES[rkey]['chain']],
            #[ATM, DNA, RESI, CH]
        #]
    #)
    #S.setRadiiWithClassifier(default_classifier)
    #SASA = sumSASA(S, freesasa.calc(S), REGEXES['DSDNA_NUCLEOTIDE_GROUPS'], 3)
    #for nkey in SASA:
        #if(nkey in NUC):
            #INT[rkey+nkey]['nuc_basa'] += (NUC[nkey]['fasa'] - SASA[nkey]['sasa'][0])/2
            #INT[rkey+nkey]['nuc_wg_basa'] += (NUC[nkey]['wg_fasa'] - SASA[nkey]['sasa'][1])/2
            #INT[rkey+nkey]['nuc_sg_basa'] += (NUC[nkey]['sg_fasa'] - SASA[nkey]['sasa'][2])/2
            #INT[rkey+nkey]['nuc_bb_basa'] += (NUC[nkey]['bb_fasa']  -SASA[nkey]['sasa'][3])/2
            #INT[rkey+nkey]['basa'] += (NUC[nkey]['fasa'] - SASA[nkey]['sasa'][0])/2

#for nkey in NUC:
    #S = buildStructure(com, 
        #in_regex=[
            #[ATM, NUC[nkey]['nucn'], NUC[nkey]['nuci'], NUC[nkey]['chain']],
            #[ATM, PRO, RESI, CH]
        #]
    #)
    #S.setRadiiWithClassifier(default_classifier)
    #SASA = sumSASA(S, freesasa.calc(S), REGEXES['RESIDUE_GROUPS'], 2)
    #for rkey in SASA:
        #if(rkey in RES):
            #INT[rkey+nkey]['res_basa'] = (RES[rkey]['fasa'] - SASA[rkey]['sasa'][0])/2
            #INT[rkey+nkey]['res_mc_basa'] = (RES[rkey]['mc_fasa'] - SASA[rkey]['sasa'][1])/2
            #INT[rkey+nkey]['res_sc_basa'] = (RES[rkey]['sc_fasa'] - SASA[rkey]['sasa'][2])/2
            #INT[rkey+nkey]['basa'] += (RES[rkey]['fasa'] - SASA[rkey]['sasa'][0])/2
