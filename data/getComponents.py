import re
import os
import json
import networkx as nx
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def compileRegexes(obj):
    # compile regexes loaded from JSON files
    objType = type(obj)
    if(objType is dict):
        for key in obj:
            if(type(obj[key]) is unicode):
                obj[key] = re.compile(obj[key])
            else:
                compileRegexes(obj[key])
    elif(objType is list):
        for i in range(len(obj)):
            if(type(obj[i]) is unicode):
                obj[i] = re.compile(obj[i])
            else:
                compileRegexes(obj[i])

def addAtomEdges(G, comp, elements):
    for i in xrange(len(comp["_chem_comp_bond.atom_id_1"])):
        atm1 = comp["_chem_comp_bond.atom_id_1"][i]
        atm2 = comp["_chem_comp_bond.atom_id_2"][i]
        if(elements[atm1] == 'H' or elements[atm2] == 'H'):
            continue
        G.add_edge(atm1, atm2)

def parseCIF(ENTRIES, cif_file):
    try:
        cif = MMCIF2Dict(cif_file)
    except:
        print('Error parsing cif file.')
        return
    
    # Skip obsolete entries
    if(cif['_chem_comp.pdbx_release_status'] == "OBS"):
        return
    
    # Only accept standard components or modified standard components
    if((
        cif['_chem_comp.type'].upper() == 'L-PEPTIDE LINKING' 
        or 
        cif['_chem_comp.type'].upper() == 'PEPTIDE LINKING'
        )
        and 
        (
        PRO.search(cif['_chem_comp.id']) 
        or 
        PRO.search(cif['_chem_comp.mon_nstd_parent_comp_id'])
        )
    ):
        hcount = cif["_chem_comp_atom.type_symbol"].count("H")
        cif["heavy_atom_count"] = len(cif["_chem_comp_atom.type_symbol"]) - hcount
        cif['_chem_comp.type'] = cif['_chem_comp.type'].upper()
        
        # Add Moieties
        mc = []
        sc = []
        for i in xrange( len(cif["_chem_comp_atom.atom_id"]) ):
            if(cif["_chem_comp_atom.type_symbol"][i] == "H"):
                # skip hydrogens
                continue
            
            if( REGEXES["PROTEIN"]["MAIN_CHAIN"].search(cif["_chem_comp_atom.atom_id"][i]) ):
                mc.append(cif["_chem_comp_atom.atom_id"][i])
            else:
                sc.append(cif["_chem_comp_atom.atom_id"][i])
        cif["main_chain_atoms"] = mc
        cif["side_chain_atoms"] = sc
        if(len(cif["main_chain_atoms"]) != 4):
            print("{}: Invalid backbone.".format(cif['_chem_comp.id']))
            return
        if(not 'CA' in cif["side_chain_atoms"]):
            print("{}: No alpha carbon found.".format(cif['_chem_comp.id']))
            return
        cif["main_chain_re"] = '|'.join(["^{}$".format(a) for a in mc])
        cif["side_chain_re"] = '|'.join(["^{}$".format(a) for a in sc])
        cif["num_main_chain"] = len(mc)
        cif["num_side_chain"] = len(sc)
        
        ENTRIES[cif['_chem_comp.id']] = cif
        print("PROTEIN: {}".format(cif['_chem_comp.id']))
    elif(cif['_chem_comp.type'].upper() == 'DNA LINKING' 
        and 
        (
        DNA.search(cif['_chem_comp.id']) 
        or
        DNA.search(cif['_chem_comp.mon_nstd_parent_comp_id'])
        )
    ):
        hcount = cif["_chem_comp_atom.type_symbol"].count("H")
        cif["heavy_atom_count"] = len(cif["_chem_comp_atom.type_symbol"]) - hcount
        cif['_chem_comp.type'] = cif['_chem_comp.type'].upper()
        
        # Get atom elements
        elements = {}
        for i in xrange(len(cif["_chem_comp_atom.type_symbol"])):
            elements[cif["_chem_comp_atom.atom_id"][i]] =  cif["_chem_comp_atom.type_symbol"][i]
        
        # Add Moieties
        G = nx.Graph()
        addAtomEdges(G, cif, elements)
        
        # Remove Sugar-non-Sugar edges
        remove = []
        
        # Find sugar-base edge
        if(cif['_chem_comp.mon_nstd_parent_comp_id'] != '?'):
            parent = cif['_chem_comp.mon_nstd_parent_comp_id']
        else:
            parent = cif['_chem_comp.id']
        
        # Determine glycosidic bond base atom
        if(REGEXES["DNA"]["PYRIMIDINES"].search(parent)):
            GB = "N1"
        elif(REGEXES["DNA"]["PURINES"].search(parent)):
            GB = "N9"
        else:
            print("{}: Unknown nucleotide type.".format(cif['_chem_comp.id']))
            return
        
        if(GB in G):
            # Remove any N1/9 - sugar bonds
            E = G.edges(GB)
            GS = None # glycosidic bond sugar atom
            for edge in E:
                s1 = (edge[0][-1] == "'") # check if this is a sugar atom
                s2 = (edge[1][-1] == "'") # check if this is a sugar atom
                if(s1):
                    GS = edge[0]
                elif(s2):
                    GS = edge[1]
            if(GS is None):
                print("{}: No glycosidic bond sugar atom found!".format(cif['_chem_comp.id']))
                return
            remove.append((GB, GS))
        else:
            print("{}: No glycosidic bond nitrogen atom found!".format(cif['_chem_comp.id']))
            return
        
        # Find sugar-phosphate edge
        Patoms = []
        for node in G.nodes():
            if(node[0] == 'P'):
                Patoms.append(node)
        if(len(Patoms) == 0):
            print("{}: No phosphate phosphorus atom found!".format(cif['_chem_comp.id']))
            return
            
        # Remove any P* - *' bonds
        for p in Patoms:
            E = G.edges(p)
            SA = None
            for edge in E:
                s1 = (edge[0][-1] == "'")
                s2 = (edge[1][-1] == "'")
                if(s1 != s2):
                    if(s1):
                        SA = edge[0]
                        remove.append((p, SA))
                    else:
                        SA = edge[1]
                        remove.append((p, SA))
        
        # Remove found edges
        for u,v in remove:
            G.remove_edge(u,v)
        
        # Use connected components to define moieties
        components = [c for c in nx.connected_components(G)]
        if(len(components) < 3):
            print("{}: Less than three moieties ({}) found!".format(cif['_chem_comp.id'], len(components)))
            return
        bs = set()
        sr = set()
        pp = set()
        unknown = []
        for c in components:
            overlap = [
                len(bs_atoms & c),
                len(sr_atoms & c),
                len(pp_atoms & c)
            ]
            maxoverlap = max(overlap)
            if(maxoverlap == 0):
                unknown.append(c)
            elif(overlap[0] == maxoverlap):
                bs = bs | c
            elif(overlap[1] == maxoverlap):
                sr = sr | c
            elif(overlap[2] == maxoverlap):
                pp = pp | c
        
        # If only one unknown, place in any empty moieties
        if(len(unknown) == 1 and len(bs) == 0 and len(sr) > 0 and len(pp) > 0):
            bs = bs | unknown[0]
            unknown = []
            print("{}: Unknown moiety ({}) added to empty base.".format(cif['_chem_comp.id'], str(unknown[0])))
        elif(len(unknown) == 1 and len(sr) == 0 and len(bs) > 0 and len(pp) > 0):
            sr = sr | unknown[0]
            unknown = []
            print("{}: Unknown moiety ({}) added to empty sugar.".format(cif['_chem_comp.id'], str(unknown[0])))
        elif(len(unknown) == 1 and len(pp) == 0 and len(bs) > 0 and len(sr) > 0):
            pp = pp | unknown[0]
            unknown = []
            print("{}: Unknown moiety ({}) added to empty phosphate.".format(cif['_chem_comp.id'], str(unknown[0])))
        
        # Classify unknown moeities
        for unk in unknown:
            stop = False
            for UA in unk:
                for rem in remove:
                    # look at removed edges
                    if(rem[0] == UA):
                        EU = UA
                        EM = rem[1]
                    elif(rem[1] == UA):
                        EU = UA
                        EM = rem[0]
                    else:
                        continue
                    # Look for edges to known moeities
                    if(EM in bs):
                        bs = bs | unk
                        stop = True
                        print("{}: Unknown moiety ({}) added to base.".format(cif['_chem_comp.id'], str(c)))
                        break
                    if(EM in sr):
                        sr = sr | unk
                        stop = True
                        print("{}: Unknown moiety ({}) added to sugar.".format(cif['_chem_comp.id'], str(c)))
                        break
                    if(EM in pp):
                        pp = pp | unk
                        stop = True
                        print("{}: Unknown moiety ({}) added to phosphate.".format(cif['_chem_comp.id'], str(c)))
                        break
                if(stop):
                    break
        
        # Check for empty moieties
        if(len(bs) == 0 or len(sr) == 0 or len(pp) == 0):
            print("{}: One or more moieties is empty!".format(cif['_chem_comp.id']))
            return
        
        # Record Sugar-Phosphate edge
        SA = None
        PA = None
        for rem in remove:
            if(rem[0] in sr and rem[1] in pp):
                SA = rem[0]
                PA = rem[1]
            elif(rem[1] in sr and rem[0] in pp):
                SA = rem[1]
                PA = rem[0]
        if(SA is None or PA is None):
            print("{}: no sugar-phosphate edge was detected!".format(cif['_chem_comp.id']))
            return
        
        # Record Sugar-Base edge
        GS = None
        GB = None
        for rem in remove:
            if(rem[0] in sr and rem[1] in bs):
                GS = rem[0]
                GB = rem[1]
            elif(rem[1] in sr and rem[0] in bs):
                GS = rem[1]
                GB = rem[0]
        if(GS is None or GB is None):
            print("{}: no sugar-base edge was detected!".format(cif['_chem_comp.id']))
            return
        
        # Write data to JSON
        print("DNA: {}".format(cif['_chem_comp.id']))
        cif["sugar-phosphate_bond"] = {
            "phosphate_atom": PA,
            "sugar_atom": SA
        }
        cif["sugar-base_bond"] = {
            "base_atom": GB,
            "sugar_atom": GS
        }
        cif["base_atoms"] = list(bs)
        cif["sugar_atoms"] = list(sr)
        cif["phosphate_atoms"] = list(pp)
        cif["base_atoms_re"] = '|'.join(["^{}$".format(a) for a in bs])
        cif["sugar_atoms_re"] = '|'.join(["^{}$".format(a) for a in sr])
        cif["phosphate_atoms_re"] = '|'.join(["^{}$".format(a) for a in pp])
        cif["num_phosphate_atoms"] = len(pp)
        cif["num_base_atoms"] = len(bs)
        cif["num_sugar_atoms"] = len(sr)
        ENTRIES[cif['_chem_comp.id']] = cif

with open('regexes.json') as FILE:
    REGEXES = json.load(FILE)
    compileRegexes(REGEXES)
DNA = REGEXES["DNA"]["STANDARD_NUCLEOTIDES"]
PRO = REGEXES["PROTEIN"]["STANDARD_RESIDUES"]

bs_atoms = set(["N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4", "O2", "O4", "C7"])
sr_atoms = set(["C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "C6'"])
pp_atoms = set(["OP1", "OP2", "OP3", "P", "P1", "P2", "PA", "PB"])

# Fix unclosed quotes
if(not os.access("components_fixed.cif", os.R_OK)):
    IN = open('components.cif')
    OUT = open('components_fixed.cif', 'w')
    qpos = [0]*50 # store up to 50 quotes
    for line in IN:
        openq = False
        qcount = 0 # count of current found quotes
        p = 0 # pointer to current character
        while(p < len(line)):
            if(line[0] == ';'):
                p += 1
                continue
            elif(line[p] == '"' or line[p] == "'"):
                quote = line[p]
                openq = True
                qpos[qcount] = p
                qcount += 1
                p += 1
                while(p < len(line)):
                    if(line[p] == quote):
                        openq = False
                        qpos[qcount] = p
                        qcount += 1
                        break
                    p += 1
            p += 1
        if(openq):
            if(qcount == 1):
                pos = qpos[0]
            else:
                pos = qpos[qcount-2]
            OUT.write(line[:pos]+line[pos+1:])
        else:
            OUT.write(line)
    IN.close()
    OUT.close()

# Break up components_fixed.cif into individual entries
ENTRIES = {}
IN = open('components_fixed.cif')
for line in IN:
    if(line[0:4] == 'data'):
        OUT = open('temp.cif', 'w')
        loop_counter = 0
        OUT.write(line)
        while(True):
            try:
                line = next(IN)
            except StopIteration:
                break
            if(line[0] == "#"):
                loop_counter += 1
                if(loop_counter > 6):
                    break
            OUT.write(line)
        OUT.close()
        parseCIF(ENTRIES, 'temp.cif')
IN.close()
os.remove("temp.cif")

OUT = open("components.json",'w')
OUT.write(json.dumps(ENTRIES,indent=None,separators=(',', ':')))
OUT.close()
os.remove("components_fixed.cif")
