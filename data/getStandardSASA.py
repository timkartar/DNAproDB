import sys
sys.path.insert(0, '..')

import freesasa
import json
import os
import numpy as np
import getBASA
from Bio.PDB.Atom import Atom
from Bio.PDB.Structure import Structure
from Bio.PDB import PDBParser
from Bio.PDB import Superimposer

from dnaprodb_utils import Regexes
from dnaprodb_utils import compileRegexes
import glob

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

def addResidue(model, rname, atoms, elements, x, y, z):
    R = model['A'][(' ', 2, ' ')]
    
    # Remove all non-backbone atoms from R
    keep = ('C', 'N', 'O', 'CA')
    delList = []
    fixList = []
    for atom in R:
        if(not atom.name in keep):
            delList.append(atom.name)
        else:
            fixList.append(atom)
    for d in delList:
        R.detach_child(d)
    
    # Change residue name
    R.resname = rname
    
    # Compute transformation 
    movList = []
    atmList = []
    for i in xrange(len(atoms)):
        if(elements[i] == 'H'):
            continue
        A = Atom(
            atoms[i],
            np.array([x[i], y[i], z[i]], dtype=np.float32),
            1.5,
            1.0,
            ' ',
            atoms[i],
            i+100,
            element=elements[i]
        )
        if(atoms[i] in keep):
            movList.append(A)
        else:
            atmList.append(A)
    sup = Superimposer()
    sup.set_atoms(fixList, movList)
    sup.apply(atmList)
    
    for i in xrange(len(atmList)):
        if(elements[i] == 'H'):
            continue
        R.add(atmList[i])

def getTripeptideSASA(structure, rname):
    OXT = 0.0
    atoms = COMPONENTS[rname]['_chem_comp_atom.atom_id']
    moieties = REGEXES.getMoietyList(rname)
    regexes = REGEXES[rname]
    SASA[rname] = {
        "total": 0.0
    }
    print(rname)
    for mty in moieties:
        SASA[rname][mty] = 0.0
    for atom in atoms:
        if(atom[0] == 'H'):
            continue
        SASA[rname][atom] = 0.0
    
    num = len(structure)
    for model in structure:
        stripH(model)
        FS = freesasa.structureFromBioPDB(model, classifier=classifier)
        sasa = freesasa.calc(FS)
        
        N = FS.nAtoms()
        for i in xrange(N):
            atom = FS.atomName(i).strip()
            if(atom == "OXT"):
                OXT += sasa.atomArea(i)/num
            resi = FS.residueNumber(i).strip()
            if(resi == "1" or resi == "3"):
                continue
            SASA[rname][atom] += sasa.atomArea(i)/num
            SASA[rname]["total"] += sasa.atomArea(i)/num
            for j in xrange(len(regexes)):
                if(regexes[j].search(atom)):
                    SASA[rname][moieties[j]] += sasa.atomArea(i)/num
                    break
    return OXT

_PROBE_RADIUS = 1.4
DATA_PATH = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(DATA_PATH,'components.json')) as FILE:
    COMPONENTS = json.load(FILE)

# Load required regexes
with open(os.path.join(DATA_PATH,'regexes.json')) as FILE:
    r = json.load(FILE)
    compileRegexes(r)
REGEXES = Regexes(regexes=r, components=COMPONENTS)

freesasa.setVerbosity(freesasa.nowarnings)
classifier = getBASA.GeneralClassifier()
classifier.initialize(DATA_PATH, COMPONENTS)
#classifier = freesasa.Classifier(fileName=os.path.join(DATA_PATH,'naccess.config'))

OXT_SASA = 0.0
COUNT = 0
parser = PDBParser(PERMISSIVE=1,QUIET=True)
SASA = {}
# Load MD simulations first
shortToLong = {
    'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS',
    'E':'GLU', 'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE',
    'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', 'P':'PRO',
    'S':'SER', 'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'
}
trajectories = glob.glob("tripeptides/*_md.pdb")
for trj in trajectories:
    rname = shortToLong[os.path.basename(trj)[1]]
    structure = parser.get_structure(rname, trj)
    OXT_SASA += getTripeptideSASA(structure, rname)
    COUNT += 1

# Now do non-standard residues
S = parser.get_structure("tripeptide", "tripeptides/tripeptide.pdb")
for component in COMPONENTS.itervalues():
    rname = component['_chem_comp.id']
    if(REGEXES["PROTEIN"]["STANDARD_RESIDUES"].search(rname)):
        continue
    if(REGEXES.isDNA(rname)):
        continue
    if(rname in SASA):
        continue
    atoms = component['_chem_comp_atom.atom_id']
    # We need a CA atom and standard backbone to proceed
    if(not ('C' in atoms and 'N' in atoms and 'O' in atoms and 'CA' in atoms)):
        print(rname)
        continue
    elements = component['_chem_comp_atom.type_symbol']
    x1 = component['_chem_comp_atom.model_Cartn_x']
    y1 = component['_chem_comp_atom.model_Cartn_y']
    z1 = component['_chem_comp_atom.model_Cartn_z']
    x2 = component['_chem_comp_atom.pdbx_model_Cartn_x_ideal']
    y2 = component['_chem_comp_atom.pdbx_model_Cartn_y_ideal']
    z2 = component['_chem_comp_atom.pdbx_model_Cartn_z_ideal']
    
    # Remove OXT - don't need this
    oi = atoms.index("OXT")
    atoms.pop(oi)
    x1.pop(oi)
    y1.pop(oi)
    z1.pop(oi)
    x2.pop(oi)
    y2.pop(oi)
    z2.pop(oi)
    elements.pop(oi)
    
    # Check for bad coordinate values
    num = 2
    skip1 = False
    skip2 = False
    if(x1.count('?') + y1.count('?') + z1.count('?') > 0):
        skip1 = True
        num -= 1
    if(x2.count('?') + y2.count('?') + z2.count('?') > 0):
        skip2 = True
        num -= 1
    if(skip1 and skip2):
        continue
    
    # Make entry in SASA dict
    SASA[rname] = {
        "total": 0.0
    }
    moieties = REGEXES.getMoietyList(rname)
    regexes = REGEXES[rname]
    for mty in moieties:
        SASA[rname][mty] = 0.0
    for i in range(len(atoms)):
        if(elements[i] == 'H'):
            continue
        SASA[rname][atoms[i]] = 0.0
    
    if(not skip1):
        addResidue(S[0], rname, atoms, elements, x1, y1, z1)
        FS = freesasa.structureFromBioPDB(S[0], classifier=classifier)
        sasa = freesasa.calc(FS, freesasa.Parameters({"probe-radius": _PROBE_RADIUS}))
        
        N = FS.nAtoms()
        for i in xrange(N):
            resi = FS.residueNumber(i).strip()
            if(resi == "1" or resi == "3"):
                continue
            atom = FS.atomName(i).strip()
            SASA[rname][atom] += sasa.atomArea(i)/num
            SASA[rname]["total"] += sasa.atomArea(i)/num
            for j in xrange(len(regexes)):
                if(regexes[j].search(atom)):
                    SASA[rname][moieties[j]] += sasa.atomArea(i)/num
                    break
    
    if(not skip2):
        addResidue(S[0], rname, atoms, elements, x2, y2, z2)
        FS = freesasa.structureFromBioPDB(S[0], classifier=classifier)
        sasa = freesasa.calc(FS, freesasa.Parameters({"probe-radius": _PROBE_RADIUS}))
        
        N = FS.nAtoms()
        for i in xrange(N):
            resi = FS.residueNumber(i).strip()
            if(resi == "1" or resi == "3"):
                continue
            atom = FS.atomName(i).strip()
            SASA[rname][atom] += sasa.atomArea(i)/num
            SASA[rname]["total"] += sasa.atomArea(i)/num
            for j in xrange(len(regexes)):
                if(regexes[j].search(atom)):
                    SASA[rname][moieties[j]] += sasa.atomArea(i)/num
                    break

# Add OXT sasa
for res in SASA:
    SASA[res]["OXT"] = OXT_SASA/COUNT

MOUT = open(os.path.join(DATA_PATH,"standard-sasa.json"),'w')
#MOUT = open(os.path.join(DATA_PATH,"standard-vasa.json"),'w')
MOUT.write(json.dumps(SASA,indent=2,separators=(',', ':'), sort_keys=True))
MOUT.close()
