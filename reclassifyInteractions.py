# NOTE: This script must be used with the same interaction distance criteria. It will not add or 
# remove any interactions, only modify moiety interactions for existing ones.
import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.3f')

import glob
import sys
import os
import re
from dnaprodb_utils import getInteractionMoiety
from dnaprodb_utils import compileRegexes
from dnaprodb_utils import Regexes
from dnaprodb_utils import C

# Get the data path
DATA_PATH = C["DATA_PATH"]

# Get interaction stats
with open(os.path.join(DATA_PATH, 'interaction_stats.json')) as FILE:
    INTERACTION_STATS = json.load(FILE)

# Load PDB Components dictionary
with open(os.path.join(DATA_PATH,'components.json')) as FILE:
    COMPONENTS = json.load(FILE)

# Load required regexes
with open(os.path.join(DATA_PATH,'regexes.json')) as FILE:
    r = json.load(FILE)
    compileRegexes(r)
REGEXES = Regexes(regexes=r, components=COMPONENTS)

def modify(json_file):
    with open(json_file) as FH:
        DATA = json.load(FH)
    
    # skip error json files
    if("error" in DATA):
        return
    
    if("num_models" not in DATA):
        return
    
    for mi in range(DATA["num_models"]):
        SSE = {}
        for sse in DATA["protein"]["models"][mi]["secondary_structure_elements"]:
            for res_id in sse["residue_ids"]:
                SSE[res_id] = sse["id"] # map res_id to parent sse_id
        
        for interface in DATA["interfaces"]["models"][mi]:
            RES_DATA = {}
            NUC_DATA = {}
            SSE_DATA = {}
            INT_STAT = {}
            for res in interface["residue_data"]:
                RES_DATA[res["res_id"]] = res
                res["interacts_by"] = set()
                res["interacts_with"] = set()
            
            for nuc in interface["nucleotide_data"]:
                NUC_DATA[nuc["nuc_id"]] = nuc
                nuc["interacts_by"] = set()
                nuc["interacts_with"] = set()
            
            for sse in interface["sse_data"]:
                SSE_DATA[sse["sse_id"]] = sse
                sse["interacts_by"] = set()
                sse["interacts_with"] = set()
            
            for stats in interface["interface_features"]:
                INT_STAT[stats["protein_chain_id"]] = stats
                stats["weak_interaction_count"] = 0
            
            # Now classify interactions and make updates
            for nr in interface["nucleotide-residue_interactions"]:
                nname = nr["nuc_name"]
                rname = nr["res_name"]
                if(not re.search(REGEXES["PROTEIN"]["STANDARD_RESIDUES"], rname)):
                    rname = COMPONENTS[rname]['_chem_comp.mon_nstd_parent_comp_id']
                if(not re.search(REGEXES["DNA"]["STANDARD_NUCLEOTIDES"], nname)):
                    nname = COMPONENTS[nname]['_chem_comp.mon_nstd_parent_comp_id']
                
                if(INTERACTION_STATS[nname+rname]["count"] >= 100):
                    # Classify interaction
                    int_mty, nuc_mty, res_mty = getInteractionMoiety(nr, INTERACTION_STATS[nname+rname], ["basa", "hbond_sum", "vdw_sum"])
                else:
                    # Try to classify interactions with fallback stats.
                    int_mty, nuc_mty, res_mty = getInteractionMoiety(nr, INTERACTION_STATS["fallback"], ["basa", "hbond_sum", "vdw_sum"])
                
                nr["moiety_interactions"] = int_mty
                nr["nucleotide_interaction_moieties"] = nuc_mty
                nr["residue_interaction_moieties"] = res_mty
                
                if(len(int_mty) == 0):
                    nr["weak_interaction"] = True
                    pro_chain = nr["res_id"][0]
                    INT_STAT[pro_chain]["weak_interaction_count"] += 1
                else:
                    nr["weak_interaction"] = False
                
                NUC_DATA[nr["nuc_id"]]["interacts_by"].update(nuc_mty)
                NUC_DATA[nr["nuc_id"]]["interacts_with"].update(res_mty)
                
                RES_DATA[nr["res_id"]]["interacts_by"].update(res_mty)
                RES_DATA[nr["res_id"]]["interacts_with"].update(nuc_mty)
                
                SSE_DATA[SSE[nr["res_id"]]]["interacts_by"].update(res_mty)
                SSE_DATA[SSE[nr["res_id"]]]["interacts_with"].update(nuc_mty)
                
            # Convert sets to lists
            for res in interface["residue_data"]:
                res["interacts_by"] = list(res["interacts_by"])
                res["interacts_with"] = list(res["interacts_with"])
            
            for nuc in interface["nucleotide_data"]:
                nuc["interacts_by"] = list(nuc["interacts_by"])
                nuc["interacts_with"] = list(nuc["interacts_with"])
            
            for sse in interface["sse_data"]:
                sse["interacts_by"] = list(sse["interacts_by"])
                sse["interacts_with"] = list(sse["interacts_with"])
    
    # Write data back to file
    with open(json_file, "w") as FH:
        FH.write(json.dumps(DATA, indent=None, separators=(',', ':'), sort_keys=True))

# loop over directories
dirlist = sys.argv[1]
DIRS = open(dirlist)
for d in DIRS:
    d = d.strip()
    for listFile in glob.glob("./{}/list*.txt".format(d)):
        LF = open(listFile)
        for pdbid in LF:
            pdbid = pdbid.strip()
            jsonFile = os.path.join(d, "{}.json".format(pdbid))
            if(os.access(jsonFile, os.R_OK)):
                modify(jsonFile)
            else:
                print(pdbid)
        LF.close()
DIRS.close()
