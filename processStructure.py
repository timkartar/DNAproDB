#!/usr/bin/env python

import numpy as np
import subprocess
import os
import re
import json
import argparse
import operator
import signal
import networkx as nx
import sys
# Biopython Disordered Atom Fix
import Bio.PDB 
copy = Bio.PDB.Atom.copy
def myCopy(self):
    shallow = copy.copy(self)
    for child in list(self.child_dict.values()):
        shallow.disordered_add(child.copy())
    return shallow
Bio.PDB.Atom.DisorderedAtom.copy=myCopy
# Biopython Disordered Atom Fix
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBParser import PDBParser 
from Bio.PDB import PDBIO
from Bio.PDB.PDBIO import Select
from Bio.Data import IUPACData
from Bio.PDB.Residue import DisorderedResidue
from Bio.PDB import NeighborSearch

# Custom Modules
import processDNA
import processProtein
import processComplex
import compileJSON
from dnaprodb_utils import log
from dnaprodb_utils import compileRegexes
from dnaprodb_utils import Regexes
from dnaprodb_utils import getIDArray
from dnaprodb_utils import getID
from dnaprodb_utils import C

from Bio import __version__ as BPV

__RESCOUNT_LIMIT = C["RESCOUNT_LIMIT"]
__TIMEOUT_LENGTH = C["TIMEOUT_LENGTH"]
__DNACOUNT_LOWER = C["DNA_MIN_COUNT"]
__PROCOUNT_LOWER = C["PRO_MIN_COUNT"]
DATA_PATH = C["DATA_PATH"]

class assembly_operation:
    def __init__(self, op_id, op_type, rotation, translation):
        self.op_id = op_id
        self.op_type = op_type
        self.rotation = np.array(rotation,dtype=np.float32).reshape(3,3)
        self.translation = np.array(translation,dtype=np.float32)

class no_hydrogen(Select):
    # Remove all hydrogens and waters
    def __init__(self, chains):
        self.chain_list = chains
    
    def accept_chain(self,chain):
        return chain.get_id() in self.chain_list
    
    #def accept_residue(self, residue):
    #    return residue.get_resname() != 'HOH'
    
    def accept_atom(self, atom):
        return (not atom.is_disordered() and atom.element != 'H')

class assembly_select(Select):
    def __init__(self, chains):
        self.chain_list = chains
    
    def accept_chain(self, chain):
        return chain.get_id() in self.chain_list
    
    def accept_atom(self, atom):
        return not atom.is_disordered()

class DNA_select(Select):
    def __init__(self, chains, components):
        self.chain_list = chains
        self.components = components
    
    def accept_chain(self, chain):
        return chain.get_id() in self.chain_list
    
    def accept_residue(self, residue):
        rname = residue.get_resname().strip()
        if(rname in self.components):
            return self.components[rname]['_chem_comp.type'] == 'DNA LINKING'
        else:
            return False
    
    def accept_atom(self, atom):
        return not atom.is_disordered()

class Protein_select(Select):
    def __init__(self, chains, components):
        self.components = components
        self.chain_list = chains
    
    def accept_chain(self, chain):
        return chain.get_id() in self.chain_list
    
    def accept_residue(self, residue):
        rname = residue.get_resname().strip()
        if(rname in self.components):
            return ( 
                self.components[rname]['_chem_comp.type'] == 'L-PEPTIDE LINKING' 
                or 
                self.components[rname]['_chem_comp.type'] == 'PEPTIDE LINKING'
            )
        else:
            return False
    
    def accept_atom(self, atom):
        return not atom.is_disordered()

def timedOut(signum, frame):
    log("Job took longer than {} seconds to complete. Aborting.".format(__TIMEOUT_LENGTH), FILE_NAME)

def assignElement(fullname):
    """Tries to guess element from atom name if not recognised."""
    name = fullname.strip()
    if name.capitalize() not in IUPACData.atom_weights:
        # Inorganic elements have their name shifted left by one position
        #  (is a convention in PDB, but not part of the standard).
        # isdigit() check on last two characters to avoid mis-assignment of
        # hydrogens atoms (GLN HE21 for example)
    
        if fullname[0].isalpha() and not (fullname[2:].isdigit() or fullname[2:] == "''"):
            putative_element = name
        else:
            # Hs may have digit in first position
            if name[0].isdigit():
                putative_element = name[1]
            else:
                putative_element = name[0]
    
        if putative_element.capitalize() in IUPACData.atom_weights:
            element = putative_element
        else:
            element = ""
        
        return element
    else:
        return name[0]

def restoreOccupancy(model, tempFile, prefix):
    # Restore occupancy, Bfactor, and element types for atoms in assembly
    ATOM_RE = re.compile('^ATOM')
    HETM_RE = re.compile('^HETATM')
    FI = open(tempFile)
    FO = open('{}_repaired.pdb'.format(prefix), 'w')
    for line in FI:
        if(ATOM_RE.search(line) or HETM_RE.search(line)):
            atm = line[12:16]
            atm_name = atm.strip()
            ch = line[21]
            resi = int(line[22:26].strip())
            ins = line[26]
            if(HETM_RE.search(line)):
                resn = line[17:20].strip()
                if(resn == 'WAT'):
                    res_id = ('W', resi, ins)
                else:
                    res_id = ('H_'+resn, resi, ins)
            else:
                res_id = (' ', resi, ins)
            
            if(model[ch].has_id(res_id) and model[ch][res_id].has_id(atm_name)):
                occupancy = model[ch][res_id][atm_name].occupancy
                bfactor = model[ch][res_id][atm_name].bfactor
                element = model[ch][res_id][atm_name].element
            else:
                occupancy = 1.0
                bfactor = 99.99
                element = assignElement(atm)
            if(re.search('O([12])P', atm_name)):
                atm_name = 'OP'+atm_name[1]
            if(len(atm_name) < 4):
                am = re.search("{}([A-Z0-9\'\*]+)?".format(element), atm_name)
                if(am.group(1)):
                    atm = "{:>2s}{:<s}".format(element, am.group(1))
                else:
                    atm = "{:>2s}  ".format(element)
            line = line[0:12]+ "{:4s}".format(atm) + line[16:54] + "{0:6.2f}{1:6.2f}          {2:>2s}{3:2s}\n".format(occupancy,bfactor,element,' ')
        
        FO.write(line)
    FI.close()
    FO.close()
    return '{}_repaired.pdb'.format(prefix)

def addAtoms(model, repaired, parser, META):
    # Read in repaired structure and add missing atoms from assembly
    repaired_assembly = parser.get_structure('repaired', repaired)
    mnum = model.get_id()
    for chain in repaired_assembly[0].get_list():
        ch = chain.get_id()
        if(ch not in model):
            continue
        for residue in chain.get_list():
            res = residue.get_id()
            if(res not in model[ch]):
                continue
            for atom in residue.get_list():
                if(atom.element == 'H'):
                    continue
                if(atom.get_id() not in model[ch][res]):
                    copy = atom.copy()
                    model[ch][res].add(copy)
                    META["added_heavy_atoms"].append("{}.{}.{}.{}".format(mnum,ch,str(res[1])+res[2].strip(),copy.get_name()))

def runPDB2PQR(pdbid, assembly, META, N):
    # Attempt to add missing heavy atoms and protonate the structure 
    # using PDB2PQR. 
    # Arguments:
    # pdbid:     Structure name/identifier
    # assembly:  Biopython structure object containing the assembly that 
    #            corresponds to pdbid.pdb
    #------------------------------------------------------------------#
    # Call PDB2PQR on the structure
    if(CLEAN_STRUCTURE):
        fileName = '{}_cleaned.pdb'.format(pdbid)
    else:
        fileName = '{}.pdb'.format(pdbid)
    
    FNULL = open(os.devnull, 'w')
    parser = PDBParser(PERMISSIVE=1,QUIET=True)
    
    if(ENSEMBLE):
        START_RE = re.compile('^MODEL')
        STOP_RE = re.compile('^ENDMDL')
        FH = open(fileName)
        i = 0
        fName = "{}_{}.pdb".format(pdbid, i)
        tName = "{}_{}.temp".format(pdbid, i)
        OUT = open(fName, "w")
        for line in FH:
            if(START_RE.search(line)):
                continue
            elif(STOP_RE.search(line)):
                OUT.close()
                rc = subprocess.call([
                        'pdb2pqr',
                        '--ff=AMBER',
                        '--keep-chain',
                        '--include-header',
                        fName,
                        tName
                    ],
                    stdout=FNULL,
                    stderr=subprocess.STDOUT
                )
                if(rc == 0 and os.access(tName, os.R_OK)):
                    repaired = restoreOccupancy(assembly[i], tName, "{}_{}".format(pdbid, i))
                    os.remove(tName)
                    addAtoms(assembly[i], repaired, parser, META)
                    os.remove(repaired)
                    os.remove(fName)
                else:
                    log("PDB2PQR failed to run, check the structure file.", pdbid)
                i += 1
                if(i == N):
                    break
                fName = "{}_{}.pdb".format(pdbid, i)
                tName = "{}_{}.temp".format(pdbid, i)
                OUT = open(fName, "w")
            else:
                OUT.write(line)
        FH.close()
    else:
        tName = '{}.temp'.format(pdbid)
        rc = subprocess.call([
                'pdb2pqr',
                '--ff=amber',
                '--chain',
                '--include-header',
                fileName,
                tName
            ],
            stdout=FNULL,
            stderr=subprocess.STDOUT
        )
        if(rc == 0 and os.access(tName, os.R_OK)):
            repaired = restoreOccupancy(assembly[0], tName, pdbid)
            os.remove(tName)
            addAtoms(assembly[0], repaired, parser, META)
            os.remove(repaired)
        else:
            log("PDB2PQR failed to run, check the structure file.", pdbid)
    
    FNULL.close()

def _disconnectGraph(G, model, threshold, removed, level="R", lengths=None):
    order = []
    if(level == 'R'):
        for node in G.nodes():
            residue = model[node[2]][node[3]]
            order.append((
                G.degree(node, weight="weight"),
                1-len(residue),
                node[2],
                node
            ))
        order.sort(key=operator.itemgetter(0,1,2), reverse=True)
    else:
        for node in G.nodes():
            order.append((
                G.degree(node, weight="weight")/float(lengths[node]),
                1-lengths[node],
                node
            ))
        order.sort(key=operator.itemgetter(0,1), reverse=True)
    
    # Greedily remove nodes until G is no longer connected
    for o in order:
        if(o[0] > threshold):
            if(level == 'R'):
                node = o[3]
                G.remove_node(node)
                residue = model[node[2]][node[3]]
                resname = residue.get_resname()
                removed[node] = {
                    "id": "{}".format(getID(residue=residue)),
                    "name": resname,
                    "reason": "many atom clashes",
                    "number_of_clashes": o[0],
                    "log": "removed clashing residue {} ({})".format(resname, getID(residue=residue))
                }
                _disconnectGraph(G, model, threshold, removed, level=level, lengths=lengths)
            else:
                node = o[2]
                G = nx.Graph(G) ## guard against frozen graph error
                G.remove_node(node)
                chain = model[node]
                for residue in chain:
                    resname = residue.get_resname()
                    removed[residue.get_full_id()] = {
                        "id": "{}".format(getID(residue=residue)),
                        "name": resname,
                        "reason": "parent chain is overlapping another chain",
                        "log": "removed residue {} ({}): parent chain is overlapping another chain".format(resname, getID(residue=residue))
                    }
                _disconnectGraph(G, model, threshold, removed, level=level, lengths=lengths)
            break

def cleanAssembly(pdbid, assembly, filter_chains, META, N):
    """Removes residues with too many missing atoms for PDB2PQR to repair,
    components which are not nucleotides or amino acids, clashing components, and fixes 
    some common atom name mistakes.
    
    Parameters
    ----------
    pdbid: string
        prefix for structure name
    assembly: Biopython structure object
        Biopython structure object containing uncleaned assembly
    filter_chains: list
        list of chain IDs to keep in the assembly
    META: dict
        META_DATA dictionary.
    N: int
        number of MODEL entries to process.
    """
    INS_CODE_LETTERS = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    with open(os.path.join(DATA_PATH,'atom-rename.json')) as FILE:
        ATOM_RENAME = json.load(FILE)

    # Remove any residue with > 50% missing atoms and fix atom names.
    # Change HETATM to ATOM for all supported residues, otherwise remove.
    REMOVED = {}
    rem_count = 0.0
    res_count = 0.0
    for i in range(N):
        for chain in assembly[i].get_list():
            chain_list = chain.get_list()
            for residue in chain.get_list():
                resname = residue.get_resname().strip()
                resid = residue.get_id()
                if(resname in COMPONENTS):
                    res_count += 1.0
                    if(resid[0][0] == 'H'):
                        # make any standard residue a non-HETATM entry
                        modified = {
                            "old_id": getID(*resid),
                            "new_id": None
                        }
                        resid = (' ', resid[1], resid[2])
                        if(resid in chain):
                            # need to modify the insertion code
                            if(resid[2] in INS_CODE_LETTERS):
                                ins = INS_CODE_LETTERS[INS_CODE_LETTERS.index(resid[2])+1]
                            else:
                                ins = 'A'
                            resid = (resid[0], resid[1], ins)
                        residue.id = resid
                        modified["new_id"] = getID(*resid)
                        META["modified_ids"].append(modified)
                    if(REGEXES.isDNA(resname)):
                        rtype = "nucleotide"
                    else:
                        rtype = "residue"
                elif(REGEXES['SOLVENT_COMPONENTS'].search(resname)):
                    # Ignore known solvent components
                    continue
                else:
                    # Delete unknown component
                    if(len(residue) == 1):
                        # except for single ions/metals
                        continue
                    
                    REMOVED[residue.get_full_id()] = {
                        "id": "{}".format(getID(residue=residue)),
                        "name": resname,
                        "reason": "unknown type",
                        "log": "removed unknown residue {} ({})".format(resname, getID(residue=residue))
                    }
                    continue
                
                if(resname in COMPONENTS):
                    if(rtype == "residue"):
                        moieties = [COMPONENTS[resname]["main_chain_re"], COMPONENTS[resname]["side_chain_re"]]
                        atom_counts = [0, 0]
                        totals = [COMPONENTS[resname]["num_main_chain"], COMPONENTS[resname]["num_side_chain"]]
                    else:
                        moieties = [COMPONENTS[resname]["base_atoms_re"], COMPONENTS[resname]["sugar_atoms_re"], COMPONENTS[resname]["phosphate_atoms_re"]]
                        atom_counts = [0, 0, 0]
                        totals = [COMPONENTS[resname]["num_base_atoms"], COMPONENTS[resname]["num_sugar_atoms"], COMPONENTS[resname]["num_phosphate_atoms"]]
                    
                    # Check for missing atoms and correct atom names
                    for atom in residue.get_list():
                        name = atom.get_fullname()
                        # Check for common misnaming of atoms
                        if(resname in ATOM_RENAME and name in ATOM_RENAME[resname]):
                            name = ATOM_RENAME[resname][name]
                            atom.fullname = name
                        # Fix atom names
                        if(len(name) < 4):
                            am = re.search("{}([A-Z0-9\'\*]+)?".format(atom.element), name)
                            if(am.group(1)):
                                name = "{:>2s}{:<s}".format(atom.element, am.group(1))
                            else:
                                name = "{:>2s}  ".format(atom.element)
                            atom.fullname = name
                        
                        # add to atom counts
                        for i in range(len(moieties)):
                            if(re.search(moieties[i], name.strip())):
                                atom_counts[i] += 1.0
                    
                    remove = False
                    if(totals[0] == 0 or totals[1] == 0):
                        REMOVED[residue.get_full_id()] = {
                            "id": "{}".format(getID(residue=residue)),
                            "name": resname,
                            "reason": "unknown type",
                            "log": "removed unknown residue {} ({})".format(resname, getID(residue=residue))
                        }
                    elif(rtype == "residue"):
                        if(atom_counts[1]/totals[1] < 0.25):
                            remove = True
                            reason = "missing too many side-chain atoms (threshold: {}%)".format(100*0.25)
                            diff = totals[1] - atom_counts[1]
                        elif(atom_counts[0]/totals[0] < 0.75):
                            remove = True
                            reason = "missing too many main-chain atoms (threshold: {}%)".format(100*0.75)
                            diff = totals[0] - atom_counts[0]
                    else:
                        if(atom_counts[0]/totals[0] < 0.75):
                            remove = True
                            reason = "missing too many base atoms (threshold: {}%)".format(100*0.75)
                            diff = totals[0] - atom_counts[0]
                        elif(atom_counts[1]/totals[1] < 0.5):
                            remove = True
                            reason = "missing too many sugar atoms (threshold: {}%)".format(100*0.5)
                            diff = totals[1] - atom_counts[1]
                    
                    # Remove incomplete residues/nucleotides
                    if(remove):
                        REMOVED[residue.get_full_id()] = {
                            "id": "{}".format(getID(residue=residue)),
                            "name": resname,
                            "reason": reason,
                            "number_missing": int(diff),
                            "log": "removed incomplete residue {} ({})".format(resname, getID(residue=residue))
                        }
            sortChain(chain)
    res_count /= N
    
    # Check for clashes and remove clashing residues
    for i in range(N):
        atom_list = []
        chain_count = {}
        for chain in assembly[i].get_list():
            cid = chain.get_id()
            chain_count[cid] = 0
            for residue in chain.get_list():
                if(REGEXES.isDNA(residue.get_resname()) or REGEXES.isProtein(residue.get_resname())):
                    chain_count[cid] += 1
                for atom in residue.get_list():
                    atom_list.append(atom)
        ns = NeighborSearch(atom_list)
        
        G = nx.Graph() # residue clashes
        clashes = ns.search_all(1.0, level='A')
        for clash in clashes:
            if(clash[0].element == "H" or clash[1].element == "H"):
                # ignore hydrogen clashes
                continue
            if(REGEXES['SOLVENT_COMPONENTS'].search(clash[0].get_parent().get_resname().strip())):
                # ignore solvent clashes
                continue
            if(REGEXES['SOLVENT_COMPONENTS'].search(clash[1].get_parent().get_resname().strip())):
                # ignore solvent clashes
                continue
            id1 = clash[0].get_parent().get_full_id()
            id2 = clash[1].get_parent().get_full_id()
            if(id1 == id2):
                # ignore self-clashes
                continue
            
            # Add edges to graph
            if(G.has_edge(id1, id2)):
                w = G.get_edge_data(id1, id2)["weight"]
                G.add_edge(id1, id2, weight=w+1)
            else:
                G.add_edge(id1, id2, weight=1)
        
        # Count chain clashes
        C = nx.Graph() # chain clashes
        clashes = ns.search_all(1.0, level='R')
        for clash in clashes:
            b1 = REGEXES.isDNA(clash[0].get_resname()) or REGEXES.isProtein(clash[0].get_resname())
            b2 = REGEXES.isDNA(clash[1].get_resname()) or REGEXES.isProtein(clash[1].get_resname())
            if(not (b1 and b2)):
                # ignore hetero stuff
                continue
            c1 = clash[0].get_parent().get_id()
            c2 = clash[1].get_parent().get_id()
            if(c1 == c2):
                # ignore self-clashes
                continue
            
            # Add edges to graph
            if(C.has_edge(c1, c2)):
                w = C.get_edge_data(c1, c2)["weight"]
                C.add_edge(c1, c2, weight=w+1)
            else:
                C.add_edge(c1, c2, weight=1)
        
        # Remove overlapping chains if any
        if(C.order() > 0):
             C_subgraphs = list(C.subgraph(c) for c in nx.connected_components(C))
             for S in C_subgraphs:
                _disconnectGraph(S, assembly[i], 0.5, REMOVED, level='C', lengths=chain_count)
        # Remove already removed nodes from G
        for r in REMOVED:
            if(G.has_node(r)):
                G.remove_node(r)
        if(G.order() > 0):
            # ARI UPDATE
            for component in nx.connected_components(G):
                S = G.subgraph(component).copy()
                _disconnectGraph(S, assembly[i], 2, REMOVED)
    
    # Verify that each MODEL has same residue/nucleotide set
    SETS = []
    for i in range(N):
        s = set()
        for chain in assembly[i]:
            for residue in chain:
                resname = residue.get_resname().strip()
                if(REGEXES.isDNA(resname) or REGEXES.isProtein(resname)):
                    s.add(getID(residue=residue))
        SETS.append(s)
    for i in range(N):
        for j in range(i+1, N):
            diff = SETS[i]^SETS[j]
            for resid in diff:
                cid, rnum, ins = resid.split('.')
                rid = (' ', int(rnum), ins)
                if(resid in SETS[i]):
                    resname = assembly[i][cid][rid].get_resname().strip()
                else:
                    resname = assembly[j][cid][rid].get_resname().strip()
                REMOVED[(assembly.get_id(), i, cid, rid)] = {
                    "id": resid,
                    "name": resname,
                    "reason": "residue appears in some models and not others",
                    "log": "removed inconsistent residue {} ({})".format(resname, resid)
                }
    
    # Check how many residues we are removing
    for r in REMOVED:
        if(REMOVED[r]["name"] in COMPONENTS):
            rem_count += 1.0
    if(rem_count/(N*res_count) > 0.2):
        log("Too many removed residues, check if something is wrong. ({}).".format(len(REMOVED)), pdbid, removed=list(REMOVED.values()))
    
    # Remove residues
    for resid in REMOVED:
        r = REMOVED[resid]
        log(r["log"], pdbid, Exit=False)
        META["removed_residues"].append(r)
        cid = resid[2]
        rid = resid[3]
        for i in range(N):
            if(rid in assembly[i][cid]):
                assembly[i][cid].detach_child(rid)
    
    POUT = PDBIO()
    POUT.set_structure(assembly)
    POUT.save("{}_cleaned.pdb".format(pdbid), assembly_select(filter_chains))

def protonate(pdbid, N):
    # Protonates the structure named pdbid-noH.pdb
    # Arguments:
    # pdbid: prefix for structure name
    #-------------------------------------------------------------------
    OUT = open("{}.pdb".format(pdbid), 'w')
    FNULL = open(os.devnull, 'w')
    rc = subprocess.call(['reduce', '-NOFLIP', '-Quiet', '-DB', os.path.join(DATA_PATH,'reduce_wwPDB_het_dict.txt'), '{}-noH.pdb'.format(pdbid)],
            stdout=OUT,
            #stderr=FNULL
        )
    FNULL.close()
    OUT.close()
    if(not os.access('{}.pdb'.format(pdbid), os.R_OK)):
            log("Reduce failed to produce any output or returned with errors.", pdbid)
    else:
        # Check that hydrogens were actually added
        FH = open("{}.pdb".format(pdbid), 'r').readlines()
        reduceRe = re.compile('^USER  MOD .+found=(\d+), std=\d+, add=(\d+)')
        for line in FH:
            rm = reduceRe.search(line)
            if(rm):
                hfound = float(rm.group(1))
                hadded = float(rm.group(2))
                if( 20*(hfound+hadded)/(len(FH)/N) > 5.0):
                    return
                else:
                    log("Structure could not be protonated.", pdbid)
        log("Structure could not be protonated.", pdbid)

def sortChain(chain):
    # determine intrinsic order of the chain
    clist = chain.get_list()
    order = 0
    for i in range(1, len(clist)):
        r0 = clist[i-1]
        r1 = clist[i]
        if(r1.get_id()[1] > r0.get_id()[1]):
            order += 1
        else:
            order -= 1
    
    # sort list
    if(order > 0):
        chain.child_list.sort()
    else:
        chain.child_list.sort(reverse=True)

def buildAssemblies(pdbid, asymmetric_unit, mmcif_dict, META):    
    # Remove any extra models from asymmetric unit if ensemble mode is 
    # turned off
    if(not ENSEMBLE):
        for i in range(1,asymmetric_unit.__len__()):
            asymmetric_unit.detach_child(i)
        N = 1
    else:
        N = len(asymmetric_unit)
    
    # Remove disordered atoms
    for i in range(N):
        for chain in asymmetric_unit[i].get_list():
            cid = chain.get_id()
            for residue in chain.get_list():
                rid = residue.get_id()
                if isinstance(residue, DisorderedResidue):
                    selected = residue.selected_child
                    selected.disordered_flag = 0
                    asymmetric_unit[i][cid].detach_child(rid)
                    asymmetric_unit[i][cid].add(selected)
                    residue = selected
                for atom in residue.get_list():
                    if atom.is_disordered():
                        selected = atom.selected_child
                        selected.altloc = " "
                        selected.disordered_flag = 0
                        asymmetric_unit[i][cid][rid].detach_child(atom.get_id())
                        asymmetric_unit[i][cid][rid].add(selected)
            # Sort chain by residue ID - appended residues may have changed in order
            sortChain(chain)
    
    if('_pdbx_struct_oper_list.id' in mmcif_dict):
        # Read the bioligcal assembly operations and apply them.
        
        # Create chain map to map chains from the asymmetric unit to those
        # of the biological assembly
        chain_map = {}
        for i in range(len(mmcif_dict['_pdbx_poly_seq_scheme.asym_id'])):
            asym_id = mmcif_dict['_pdbx_poly_seq_scheme.asym_id'][i]
            chain_id = mmcif_dict['_pdbx_poly_seq_scheme.pdb_strand_id'][i]
            if(not asym_id in chain_map):
                chain_map[asym_id] = chain_id
        if('_pdbx_nonpoly_scheme.asym_id' in mmcif_dict):
            for i in range(len(mmcif_dict['_pdbx_nonpoly_scheme.asym_id'])):
                asym_id = mmcif_dict['_pdbx_nonpoly_scheme.asym_id'][i]
                chain_id = mmcif_dict['_pdbx_nonpoly_scheme.pdb_strand_id'][i]
                if(not asym_id in chain_map):
                    chain_map[asym_id] = chain_id
        
        # Extract Biological Assembly Operations
        operation_ids = mmcif_dict['_pdbx_struct_oper_list.id']
        assembly_operations = {} # dictionary which stores assembly_operation objects
        if(type(operation_ids) is list):
            # Multiple biological assembly operations
            for i in range(len(operation_ids)):
                assembly_operations[operation_ids[i]] = assembly_operation(
                    operation_ids[i],
                    mmcif_dict['_pdbx_struct_oper_list.type'][i],
                    [
                        mmcif_dict['_pdbx_struct_oper_list.matrix[1][1]'][i],
                        mmcif_dict['_pdbx_struct_oper_list.matrix[1][2]'][i],
                        mmcif_dict['_pdbx_struct_oper_list.matrix[1][3]'][i],
                        mmcif_dict['_pdbx_struct_oper_list.matrix[2][1]'][i],
                        mmcif_dict['_pdbx_struct_oper_list.matrix[2][2]'][i],
                        mmcif_dict['_pdbx_struct_oper_list.matrix[2][3]'][i],
                        mmcif_dict['_pdbx_struct_oper_list.matrix[3][1]'][i],
                        mmcif_dict['_pdbx_struct_oper_list.matrix[3][2]'][i],
                        mmcif_dict['_pdbx_struct_oper_list.matrix[3][3]'][i],
                    ],
                    [
                        mmcif_dict['_pdbx_struct_oper_list.vector[1]'][i],
                        mmcif_dict['_pdbx_struct_oper_list.vector[2]'][i],
                        mmcif_dict['_pdbx_struct_oper_list.vector[3]'][i],
                        
                    ],
                )
        else:
            assembly_operations[operation_ids] = assembly_operation(
                operation_ids,
                mmcif_dict['_pdbx_struct_oper_list.type'],
                [
                    mmcif_dict['_pdbx_struct_oper_list.matrix[1][1]'],
                    mmcif_dict['_pdbx_struct_oper_list.matrix[1][2]'],
                    mmcif_dict['_pdbx_struct_oper_list.matrix[1][3]'],
                    mmcif_dict['_pdbx_struct_oper_list.matrix[2][1]'],
                    mmcif_dict['_pdbx_struct_oper_list.matrix[2][2]'],
                    mmcif_dict['_pdbx_struct_oper_list.matrix[2][3]'],
                    mmcif_dict['_pdbx_struct_oper_list.matrix[3][1]'],
                    mmcif_dict['_pdbx_struct_oper_list.matrix[3][2]'],
                    mmcif_dict['_pdbx_struct_oper_list.matrix[3][3]'],
                ],
                [
                    mmcif_dict['_pdbx_struct_oper_list.vector[1]'],
                    mmcif_dict['_pdbx_struct_oper_list.vector[2]'],
                    mmcif_dict['_pdbx_struct_oper_list.vector[3]'],
                    
                ],
            )
        
        # Build Biological Assemblies
        assembly_ids = mmcif_dict['_pdbx_struct_assembly_gen.assembly_id']
        if(not type(assembly_ids) is list):
            assembly_ids = [assembly_ids]
            mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list'] = [mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list']]
            mmcif_dict['_pdbx_struct_assembly_gen.oper_expression'] = [mmcif_dict['_pdbx_struct_assembly_gen.oper_expression']]
        
        # Check for ranges in .oper_expression list
        for i in range(len(mmcif_dict['_pdbx_struct_assembly_gen.oper_expression'])):
            op = mmcif_dict['_pdbx_struct_assembly_gen.oper_expression'][i]
            m = re.search("(\d+)-(\d+)", op)
            if(m is not None):
                s = int(m.group(1))
                e = int(m.group(2)) + 1
                op = ",".join([str(_) for _ in range(s, e)])
                mmcif_dict['_pdbx_struct_assembly_gen.oper_expression'][i] = op
        
        # create intructions for each assembly (keyed by .assembly_id)
        assembly_instructions = {}
        for i in range(len(assembly_ids)):
            if(assembly_ids[i] in assembly_instructions):
                assembly_instructions[assembly_ids[i]]['operations'] += mmcif_dict['_pdbx_struct_assembly_gen.oper_expression'][i].split(',')
                for j in range(len(mmcif_dict['_pdbx_struct_assembly_gen.oper_expression'][i].split(','))):
                    assembly_instructions[assembly_ids[i]]['op_chains'].append(mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list'][i].split(','))
            else:
                assembly_instructions[assembly_ids[i]] = {
                    'operations': mmcif_dict['_pdbx_struct_assembly_gen.oper_expression'][i].split(','),
                    'op_chains': []
                }
                for j in range(len(mmcif_dict['_pdbx_struct_assembly_gen.oper_expression'][i].split(','))):
                    assembly_instructions[assembly_ids[i]]['op_chains'].append(mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list'][i].split(','))
        
        # convert .asym_id_list to chain ids
        for instruction in assembly_instructions.values():
            for i in range(len(instruction['op_chains'])):
                chain_ids = []
                for asym_id in instruction['op_chains'][i]: #guard against chain_map key error
                    if(asym_id in chain_map.keys() and (not chain_map[asym_id] in chain_ids)):
                        chain_ids.append(chain_map[asym_id])
                instruction['op_chains'][i] = chain_ids
        
        # Generate each assembly
        for key in sorted(assembly_instructions.keys()):
            id_list = [
                'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
                'P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d',
                'e','f','g','h','i','j','k','l','m','n','o','p','q','r','s',
                't','u','v','w','x','y','z','0','1','2','3','4','5','6','7',
                '8','9'
            ]
            
            # compute number of generated chains
            chain_count = 0
            for i in range(len(assembly_instructions[key]['operations'])):
                chain_count += len(assembly_instructions[key]['op_chains'][i])
            if(chain_count > len(id_list)):
                log("This structure is too large, aborting. ({} chains)".format(chain_count), pdbid, chain_count=chain_count)
            
            # create a blank structure
            assembly = asymmetric_unit.copy()
            for i in range(N):
                chain_ids = []
                for chain in assembly[i].child_list:
                    chain_ids.append(chain.get_id())
                for chain_id in chain_ids:
                    assembly[i].detach_child(chain_id)
        
            # apply operations to each chain
            filter_chains = []
            for i in range(len(assembly_instructions[key]['operations'])):
                op_id = assembly_instructions[key]['operations'][i]
                for chain_id in assembly_instructions[key]['op_chains'][i]:
                    for j in range(N):
                        if(not asymmetric_unit[j].has_id(chain_id)):
                            continue
                        chain = asymmetric_unit[j][chain_id].copy()
                        chain.transform(
                            assembly_operations[op_id].rotation, 
                            assembly_operations[op_id].translation
                        )
                        # rename all chains to be a single letter
                        if(len(chain.get_id()) > 1):
                            for c in id_list:
                                try:
                                    chain.id = c 
                                    break
                                except(ValueError):
                                    pass
                        # check if assembly already contains this chain
                        if(assembly[j].has_id(chain.get_id())):
                            for cid in id_list:
                                if(not assembly[j].has_id(cid)):
                                    chain.id = cid 
                                    break
                        assembly[j].add(chain)
                    filter_chains.append(chain.get_id())
                    if(chain.get_id() in id_list):
                        id_list.remove(chain.get_id())
                    META["assembly_chains"][chain.get_id()] = chain_id
            break # currently only process first assembly
    else:
        # This structure lacks any biological assembly instructions, 
        # most likely an NMR structure.
        filter_chains = []
        assembly =  asymmetric_unit
        for i in range(N):
            for chain in assembly[i].child_list:
                cid = chain.get_id()
                if(cid not in filter_chains):
                    filter_chains.append(cid)
                META["assembly_chains"][cid] = cid
    
    return assembly, filter_chains, N

def processPDBFile(pdbid, pdb_file, META):
    parser = PDBParser(PERMISSIVE=1,QUIET=True)
    assembly = parser.get_structure(pdbid, pdb_file)
    chain_ids = []
    if(ENSEMBLE):
        N = len(assembly)
    else:
        N = 1
    
    # Remove disordered atoms
    for i in range(N):
        for chain in assembly[i].get_list():
            cid = chain.get_id()
            for residue in chain.get_list():
                rid = residue.get_id()
                if isinstance(residue, DisorderedResidue):
                    selected = residue.selected_child
                    selected.disordered_flag = 0
                    assembly[i][cid].detach_child(rid)
                    assembly[i][cid].add(selected)
                    residue = selected
                for atom in residue.get_list():
                    if atom.is_disordered():
                        selected = atom.selected_child
                        selected.altloc = " "
                        selected.disordered_flag = 0
                        assembly[i][cid][rid].detach_child(atom.get_id())
                        assembly[i][cid][rid].add(selected)
            # Sort chain - appended residues may have changed in order
            sortChain(chain)
    
    for i in range(N):
        for chain in assembly[i].child_list:
            cid = chain.get_id()
            if(cid not in chain_ids):
                chain_ids.append(cid)
            META["assembly_chains"][cid] = cid
    
    return assembly, chain_ids, N

def writeStructures(pdbid, assembly, filter_chains, COMPONENTS, c=True, d=True, p=True):
    io = PDBIO()
    io.set_structure(assembly)
    if(c):
        cname = "{}.pdb".format(pdbid)
        io.save(cname, assembly_select(filter_chains))
    if(d):
        dname = "{}-DNA.pdb".format(pdbid)
        io.save(dname, DNA_select(filter_chains, COMPONENTS))
    if(p):
        pname = "{}-protein.pdb".format(pdbid)
        io.save(pname, Protein_select(filter_chains, COMPONENTS))

def main(file_name, mPRE_PDB2PQR=False):
    """This module processes a PDB or mmCIF file and performs the 
    following functions:
    1. Builds the biological assembly for a given mmCIF file. If a PDB 
        file is given, then it is assumed to already contain the 
        biological assembly coordinates.
    2. 
    3.
    4.
    
    Parameters
    ----------
    pdbid: string
        A PDB id or prefex of the file to be processed. Should be named 
        either 'pdbid'.pdb or 'pdbid'.cif. 
    """
    PRE_PDB2PQR = mPRE_PDB2PQR
    signal.signal(signal.SIGALRM, timedOut)
    signal.alarm(__TIMEOUT_LENGTH)
    pdbid, ext = file_name.split('.')
    print(("Structure ID: {}".format(pdbid)))
    print(("BioPython Version: {}".format(BPV)))
    cif_file = "{}.cif".format(pdbid)
    pdb_file = "{}.pdb".format(pdbid)
    
    # Store various data needed later
    META_DATA = {
        "assembly_chains": {},
        "removed_residues": [],
        "added_heavy_atoms":[],
        "modified_ids": [],
        "deleted_models": [],
        "options":{
            "run_pdb2pqr": PRE_PDB2PQR,
            "clean_structure": CLEAN_STRUCTURE,
            "ensemble_analysis": ENSEMBLE,
            "include_annotations": ADD_MMCIF
        },
        "program_versions": {
            "dssr": "1.6.9",
            "curves": "5.3",
            "x3dna": "2.3",
            "dssp": "2.0.4",
            "hbplus": "3.2",
            "reduce": "3.24.130724",
            "pdb2pqr": "2.1.1",
            "msms": "2.6.1",
            "x3dna-snap": "beta-r10-2017apr10"
        }
    }
    if ext == "cif":
        # Process mmCIF file if exists
        parser = MMCIFParser(QUIET=True)
        asymmetric_unit = parser.get_structure(pdbid, file_name)
        mmcif_dict = MMCIF2Dict(file_name)
        
        # Create the biological assembly(s)
        assembly, filter_chains, N = buildAssemblies(pdbid, asymmetric_unit, mmcif_dict, META_DATA)
        PDB = False
    elif(os.access(pdb_file, os.R_OK)):
        # Process PDB file if exists
        assembly, filter_chains, N = processPDBFile(pdbid, file_name, META_DATA)
        PDB = True
        mmcif_dict = None
    else:
        print("No PDB/mmCIF file found.", pdbid)
        return "No PDB/mmCIF file found."
    META_DATA["ensemble"] = ENSEMBLE and (N > 1)
    
    # Clean the assembly
    if(CLEAN_STRUCTURE):
        print("YO I HERE")
        cleanAssembly(pdbid, assembly, filter_chains, META_DATA, N)
    elif(not PDB):
        io = PDBIO()
        io.set_structure(assembly)
        io.save("{}.pdb".format(pdbid), assembly_select(filter_chains))
    
    # Check size of the assembly
    dnacount = 0
    procount = 0
    for chain in assembly[0]:
        for residue in chain:
            resname = residue.get_resname()
            if(REGEXES.isDNA(resname)):
                dnacount += 1
            if(REGEXES.isProtein(resname)):
                procount += 1
    if((dnacount+procount) > __RESCOUNT_LIMIT):
        print("This structure is too large, aborting. ({} components)".format(procount+dnacount), pdbid, component_count=dnacount+procount)
        return
    elif(dnacount < __DNACOUNT_LOWER):
        print("This structure contains too few nucleotides. ({} nucleotides)".format(dnacount), pdbid, nuc_count=dnacount)
        return
    elif(procount < __PROCOUNT_LOWER):
        print("This structure contains too few residues. ({} residues)".format(procount), pdbid, res_count=procount)
        return
     
    # Optionally repair with PDB2PQR
    if(PRE_PDB2PQR):
        print("RUNNING PDB2PQR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        runPDB2PQR(pdbid, assembly, META_DATA, N)
    else:
        print("NOOOT RUNNING PDB2PQR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("Past running PDB2PQR")

    if(os.access('{}_cleaned.pdb'.format(pdbid), os.R_OK)):
        os.remove('{}_cleaned.pdb'.format(pdbid))
    print("Here 1")
    if(not ADD_MMCIF):
        mmcif_dict = None
    
    # Get list of residue/nucleotide IDs
    IDs = getIDArray(assembly[0], REGEXES)
    if(len(IDs["protein"]) == 0):
        print("No protein residues found in biological assembly!", pdbid)
        return
    elif(len(IDs["protein"]) < C["MINIMUM_RES_COUNT"]):
        print("Less than {} protein residues found - aborting.".format(C["MINIMUM_RES_COUNT"]), pdbid)
        return
    elif(len(IDs["dna"]) == 0):
        print("No DNA nucleotides found in biological assembly!", pdbid)
        return
    print("Here 2")
    # Write out DNA 
    writeStructures(pdbid, assembly, filter_chains, COMPONENTS, c=False, p=False)
    print("Here 3")
    # Process DNA
    DNA_DATA, REMOVE = processDNA.process("{}-DNA".format(pdbid), N, REGEXES, COMPONENTS, META_DATA)
    print("Here 4") 
    if(len(REMOVE) > 0):
        # remove nucleotides not recognized by DSSR
        for nid in REMOVE:
            cid, num, ins = nid.split('.')
            rid = (' ', int(num), ins)
            for i in range(N):
                resname = assembly[i][cid][rid].get_resname()
                assembly[i][cid].detach_child(rid)
                log("removed unrecognized nucleotide {} ({})".format(resname, nid), pdbid, Exit=False, name=resname)
                META_DATA["removed_residues"].append({
                    "id": "{}".format(nid),
                    "name": resname,
                    "reason": "not recognized by dssr",
                    "model": i
                })
    
    NUCLEOTIDES = []
    DELETE_MODELS = []
    for i in range(len(DNA_DATA)):
        # Check for empty entities
        if(len(DNA_DATA[i]["entities"]) == 0):
            DELETE_MODELS.append(i)
            N -= 1
            continue
        
        # Loop over each model entry
        NUCLEOTIDES.append({})
        for n in DNA_DATA[i]["nucleotides"]:
            # compile groove atoms regexes
            if("groove_atoms" in n):
                n["groove_atoms"]["sg"] = re.compile(
                    '|'.join(["^{}$".format(_) for _ in n["groove_atoms"]["sg"]])
                )
                n["groove_atoms"]["wg"] = re.compile(
                    '|'.join(["^{}$".format(_) for _ in n["groove_atoms"]["wg"]])
                )
            NUCLEOTIDES[-1][n["id"]] = n
    
    if(N == 0):
        # this structure contains no DNA entities - throw error and exit
        log("No valid DNA entities were found!", pdbid)
    
    if(len(DELETE_MODELS) > 0):
        # delete models with no DNA entities
        for i in DELETE_MODELS:
            DNA_DATA.pop(i)
            assembly.detach_child(i)
            META_DATA["deleted_models"].append({
                "model": i,
                "reason": "No suitable DNA entities were detected."
            })
            log("deleted model {}, no DNA entities detected.".format(i), pdbid, Exit=False)
    REGEXES.nucleotides = NUCLEOTIDES
    if(DEBUG_DNA):
        exit(0)
    
    # Write unprotonated structure to file
    io = PDBIO()
    io.set_structure(assembly)
    hname = "{}-noH.pdb".format(pdbid)
    io.save(hname, no_hydrogen(filter_chains))
    
    # Protonate Structure with Reduce
    protonate(pdbid, N)
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    assembly = parser.get_structure(pdbid, "{}.pdb".format(pdbid))
    
    # Write protonated complex, protein and DNA to file (models may have been removed from original assembly)
    writeStructures(pdbid, assembly, filter_chains, COMPONENTS)
    
    # Process Protein
    PRO_DATA, DSSP = processProtein.process(pdbid, N, COMPONENTS, REGEXES, IDs,
        meta_data=META_DATA, sse_type='consensus', mmcif_dict=mmcif_dict
    )
    
    # Process Interactions
    INT_DATA = processComplex.process(pdbid, N, COMPONENTS, assembly, DSSP, DATA_PATH, REGEXES, NUCLEOTIDES, IDs)
    
    # Combine all data for final JSON output
    compileJSON.comp(pdbid, N, assembly, PRO_DATA, DSSP, DNA_DATA, INT_DATA, REGEXES, NUCLEOTIDES, COMPONENTS, META_DATA,
        mmcif_dict=mmcif_dict, ADD_MMCIF=ADD_MMCIF
    )
    
    # Write out the meta-data
    MOUT = open("{}-meta.json".format(pdbid),'w')
    MOUT.write(json.dumps(META_DATA,indent=None,separators=(',', ':')))
    MOUT.close()
    signal.alarm(0)
    return True


########################################################################
# Get arguments
# parser = argparse.ArgumentParser()
# parser.add_argument("file_name",
#                 help="PDB or mmCIF file of a protein-DNA complex structure.")
# parser.add_argument("-p", "--pdb2pqr", dest="pdb2pqr", action='store_true',
#                 help="Process the structure with PDB2PQR.")

# parser.add_argument("-D", "--DNA_debug", dest="debug_dna", action='store_true',
#                 help="Exit after processing DNA (debug).")

# group1 = parser.add_mutually_exclusive_group()
# group1.add_argument("-c", "--clean", dest="clean", action='store_true',
#                 help="If -c|--clean is selected, attempt to protonate and add missing heavy atoms.")
# group1.add_argument("-C", "--no_clean", dest="clean", action='store_false',
#                 help="If -C|--no_clean is selected, do not attempt to protonate and add missing heavy atoms.")

# group2 = parser.add_mutually_exclusive_group()
# group2.add_argument("-e", "--ensemble", dest="ensemble", action='store_true',
#                 help="Treat multiple MODEL entries as a statistical ensemble, if present.")
# group2.add_argument("-E", "--no_ensemble", dest="ensemble", action='store_false',
#                 help="Ignore multiple MODEL entries, only process the first.")

# group3 = parser.add_mutually_exclusive_group()
# group3.add_argument("-m", "--meta", dest="meta", action='store_true',
#                 help="Process meta data from mmcif file.")
# group3.add_argument("-M", "--no_meta", dest="meta", action='store_false',
#                 help="Do not attempt to process meta data from mmcif file.")

# parser.set_defaults(clean=True, pdb2pqr=False, ensemble=True, meta=False, debug_dna=False)
# args = parser.parse_args()

# FILE_NAME = args.file_name
# CLEAN_STRUCTURE = args.clean
# PRE_PDB2PQR = args.pdb2pqr
# ENSEMBLE = args.ensemble
# ADD_MMCIF = args.meta
# DEBUG_DNA = args.debug_dna


FILE_NAME = ""
CLEAN_STRUCTURE = True
PRE_PDB2PQR = False
ENSEMBLE = True
ADD_MMCIF = False
DEBUG_DNA = False

# Load PDB Components dictionary
with open(os.path.join(DATA_PATH,'components.json')) as FILE:
    COMPONENTS = json.load(FILE)

# Load required regexes
with open(os.path.join(DATA_PATH,'regexes.json')) as FILE:
    r = json.load(FILE)
    compileRegexes(r)
REGEXES = Regexes(regexes=r, components=COMPONENTS)

if __name__ == '__main__':
    main(sys.argv[1])
