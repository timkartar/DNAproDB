import re
import numpy as np
import json
import freesasa
from Bio.PDB import PDBIO
import os

# Identifier field regexes
CHAIN_RE = '[a-zA-Z0-9]'
RESI_RE = '-?\d+'
RESN_RE = '[A-Z]+'
ATOM_RE = '[A-Z]+[0-9]*[\'\*]*'
INS_RE = '[A-Z\s]'
ID_RE = "{}\.{}\.{}".format(CHAIN_RE, RESI_RE, INS_RE)


# Load configuration file
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(BASE_DIR, "config.json")) as FH:
    C = json.load(FH)
C["BASE_DIR"] = BASE_DIR

class Regexes(object):
    def __init__(self, regexes=None, components=None, nucleotides=None, i=0):
        self.components = components
        self.nucleotides = nucleotides
        self.regexes = regexes
        self.i = i # Default model number
    
    def __getitem__(self, key):
        """ Return a regex or list or regexes based on the value of 
        'key'.
        If key is in self.regexes, return a compiled regex from the list
        of pre-defined regexes.
        If key is in self.nucleotides, return a list of regexes based on
        the nucleotide secondary structure corresponding to the 
        nucleotide moeities.
        If key is in self.components, assume this is a peptide-linking 
        residue and return a list of regexes corresponding the the main 
        chain and side chain.
        """ 
        
        if(self.regexes is not None and key in self.regexes):
            # return a pre-defined regex or object in regexes dictionary
            return self.regexes[key]
        elif(self.nucleotides is not None and key in self.nucleotides[self.i]):
            # return a list of regexes based on nucleotide secondary 
            # structure assignment
            ss = self.nucleotides[self.i][key]["secondary_structure"]
            name = self.nucleotides[self.i][key]["name"]
            
            # Compile regexes if not already so
            if(not isinstance(self.components[name]["base_atoms_re"], re.Pattern)):
                self.components[name]["base_atoms_re"] = re.compile(self.components[name]["base_atoms_re"])
            if(not isinstance(self.components[name]["sugar_atoms_re"], re.Pattern)):
                self.components[name]["sugar_atoms_re"] = re.compile(self.components[name]["sugar_atoms_re"])
            if(not isinstance(self.components[name]["phosphate_atoms_re"], re.Pattern)):
                self.components[name]["phosphate_atoms_re"] = re.compile(self.components[name]["phosphate_atoms_re"])
                
            if(ss == "helical"):
                return [
                    self.nucleotides[self.i][key]["groove_atoms"]["wg"],
                    self.nucleotides[self.i][key]["groove_atoms"]["sg"],
                    self.components[name]["sugar_atoms_re"],
                    self.components[name]["phosphate_atoms_re"]
                ]
            else:
                return [
                    self.components[name]["base_atoms_re"],
                    self.components[name]["sugar_atoms_re"],
                    self.components[name]["phosphate_atoms_re"]
                ]
        elif(self.components is not None and key in self.components):
            # return backbone and side-chain regexes
            return [
                self.regexes["PROTEIN"]["MAIN_CHAIN"],
                self.regexes["PROTEIN"]["SIDE_CHAIN"]
            ]
        else:
            raise ValueError("key: '{}' not recognized!".format(key))
    
    def __contains__(self, key):
        if((self.regexes is not None) and (key in self.regexes)):
            return True
        elif((self.nucleotides is not None) and (key in self.nucleotides[self.i])):
            return True
        elif((self.components is not None) and (key in self.components)):
            return True
        return False
    
    def isProtein(self, resname):
        resname = resname.strip()
        if(resname in self.components):
            return (
                self.components[resname]['_chem_comp.type'] == 'L-PEPTIDE LINKING' 
                or 
                self.components[resname]['_chem_comp.type'] == 'PEPTIDE LINKING'
            )
        return False
    
    def isDNA(self, resname):
        resname = resname.strip()
        if(resname in self.components):
            return self.components[resname]['_chem_comp.type'] == 'DNA LINKING'
        return False
    
    def getMoietyList(self, key):
        if(self.nucleotides is not None and key in self.nucleotides[self.i]):
            # return a list of regexes based on nucleotide secondary 
            # structure assignment
            ss = self.nucleotides[self.i][key]["secondary_structure"]
            name = self.nucleotides[self.i][key]["name"]
            
            if(ss == "helical"):
                return C["NUC_MTY_LABEL_DS"]
            else:
                return C["NUC_MTY_LABEL_SS"]
        
        elif(self.components is not None and key in self.components):
            # return backbone and side-chain regexes
            return C["RES_MTY_LABEL"]
    
    def setModel(self, i):
        self.i = i

def log(message, pdbid, Exit=True, print_json=True, trace=False, **kwargs):
    """Crude error logging function.
    
    Parameters
    ----------
    message: string
        The message to be logged.
    pdbid: string
        A pdbid or prefix to indicate which file/structure is creating
        the issue.
    Exit: boolean
        Wether to halt program execution or not. Warnings will generally
        use Exit=False and errors Exit=True.
    """
    LOG = open('logfile.txt', 'a')
    LOG.write('{}: {}\n'.format(pdbid, message))
    LOG.close()
    
    if(Exit):
        if(print_json):
            JSON = {
                "message": message,
                "pdbid": pdbid,
                "error": "processing error"
            }
            for kw in kwargs:
                JSON[kw] = kwargs[kw]
            EOUT = open("{}.json".format(pdbid), "w")
            EOUT.write(json.dumps(JSON, indent=None, separators=(',', ':')))
            EOUT.close()
        exit(1)

def getInteractionMoiety(interaction, stats, fields):
    moiety_labels = C["NUC_MTY_LABEL"] + C["RES_MTY_LABEL"]
    # Compute ratios based on statistics
    ratios = {}
    for f in fields:
        for mty1 in interaction[f]:
            if(mty1 not in moiety_labels):
                continue
            if(mty1 not in ratios):
                # add dict for mty1
                ratios[mty1] = {}
            for mty2 in interaction[f][mty1]:
                if(mty2 not in moiety_labels):
                    continue
                if(mty2 not in ratios[mty1]):
                    # add accumulator for mty1-mty2
                    ratios[mty1][mty2] = []
                lower = stats[f][mty1][mty2]["p20"]
                upper = stats[f][mty1][mty2]["p80"]
                if(upper == 0.0):
                    ratios[mty1][mty2].append(min(1.0, interaction[f][mty1][mty2] - 1e-5))
                    continue
                elif(lower == upper):
                    denom = upper
                    shift = 1e-5
                else:
                    denom = upper - lower
                    shift = lower
                ratios[mty1][mty2].append((interaction[f][mty1][mty2] - shift)/(denom + 1e-5))
    
    # Determine interaction moieties
    moieties = set()
    for mty1 in ratios:
        for mty2 in ratios[mty1]:
            if(len(ratios[mty1][mty2]) == 0):
                continue
            if(max(ratios[mty1][mty2]) >= 0.0):
                moieties.add("{}.{}".format(mty1, mty2))
        
    # Split up moieties
    nuc_moieties = set()
    res_moieties = set()
    for mty in moieties:
        m1, m2 = mty.split('.')
        if(m2 in C["RES_MTY_LABEL"]):
            res_moieties.add(m2)
        if(m1 in C["NUC_MTY_LABEL"]):
            nuc_moieties.add(m1)
    
    return list(moieties), list(nuc_moieties), list(res_moieties)

def compileRegexes(obj):
    # compile regexes loaded from JSON files
    objType = type(obj)
    if(objType is dict):
        for key in obj:
            if(type(obj[key]) is str):
                obj[key] = re.compile(obj[key])
            else:
                compileRegexes(obj[key])
    elif(objType is list):
        for i in range(len(obj)):
            if(type(obj[i]) is str):
                obj[i] = re.compile(obj[i])
            else:
                compileRegexes(obj[i])

def getCM(residue):
    CM = np.zeros(3)
    count = 0.0
    for atom in residue:
        if(atom.element == 'H'):
            continue
        CM += atom.get_coord()
        count += 1.0
    return CM/count

def getStructureFromModel(model, classifier=None):
    outFile = "gsfm.temp.pdb"
    io = PDBIO()
    io.set_structure(model)
    io.save(outFile)
    
    if(classifier is not None):
        structure = freesasa.Structure(outFile, classifier=classifier)
    else:
        structure = freesasa.Structure(outFile)
    
    if(os.access(outFile, os.R_OK)):
        os.remove(outFile)
    
    return structure

def getIDArray(model, REGEXES):
    COMPLEX_IDs = {}
    DNA_IDs = {}
    PROTEIN_IDs = {}
    for chain in model:
        for residue in chain:
            name = residue.get_resname()
            rid = getID(residue=residue)
            ch, num, ins = rid.split('.')
            alt_id = getID(ch, num.strip(), name.strip())
            if(REGEXES.isDNA(name)):
                DNA_IDs[rid] = alt_id
                COMPLEX_IDs[rid] = alt_id
            elif(REGEXES.isProtein(name)):
                PROTEIN_IDs[rid] = alt_id
                COMPLEX_IDs[rid] = alt_id
    return {'complex': COMPLEX_IDs, 'dna': DNA_IDs, 'protein': PROTEIN_IDs}

def getID(*args, **kwargs):
    if('residue' in kwargs):
        ch = kwargs['residue'].get_parent().get_id()
        rid = kwargs['residue'].get_id()
        fields = [ch, str(rid[1]), rid[2]]
        return '.'.join(fields)
    else:
        args = [str(_) for _ in args]
        return '.'.join(args)

def roundFloats(dictionary, precision=3):
    for key, value in dictionary.items():
        if isinstance(value, float):
            dictionary[key] = round(value, precision)
        elif isinstance(value, np.float32):
            dictionary[key] = round(float(value), precision)
        elif isinstance(value, dict):
            roundFloats(value,precision)
        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    roundFloats(item,precision)

def nucleotideMoiety(atom, nuc_id, REGEXES):
    classes = REGEXES[nuc_id]
    moieties = REGEXES.getMoietyList(nuc_id)
    for i in range(len(classes)):
        if(classes[i].search(atom)):
            return moieties[i]
            break
    # atom not recognized - log error and exit
    pass

def residueMoiety(atom, resn, REGEXES):
    classes = REGEXES[resn]
    moieties = REGEXES.getMoietyList(resn)
    for i in range(len(classes)):
        if(classes[i].search(atom)):
            return moieties[i]
            break
    # atom not recognized - log error and exit
    pass

def getHash(*args):
    args = list(args)
    # Check that no tuples are present
    for item in args:
        if(isinstance(item, tuple)):
            for t in item:
                args.append(t)
            args.remove(item)
    return '@'.join(sorted(args))
