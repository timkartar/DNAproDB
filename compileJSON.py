# Python Modules
import re
import numpy as np
import math as ma
import os
import json
import datetime
import subprocess
import copy

# Custom modules
from dnaprodb_utils import getHash, log, getID, getInteractionMoiety
from dnaprodb_utils import C

from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.3f')

# Get the data path
DATA_PATH = C["DATA_PATH"]

# Get standard SASA values
with open(os.path.join(DATA_PATH,'standard-sasa.json')) as FILE:
    STANDARD_SASA = json.load(FILE)

# Get interaction stats
with open(os.path.join(DATA_PATH, 'interaction_stats.json')) as FILE:
    INTERACTION_STATS = json.load(FILE)

# Get moiety labels
nucDSMtyLabel = C["NUC_MTY_LABEL_DS"]
nucSSMtyLabel = C["NUC_MTY_LABEL_SS"]
nucMtyLabel = C["NUC_MTY_LABEL"]
resMtyLabel = C["RES_MTY_LABEL"]
resSST = C["RES_SST"]

# Get CV cutoffs
CV_PEAK_CUTOFF = lambda x: 0.00 <= x < C["CV_PEAK_UPPER"]
CV_FLAT_CUTOFF = lambda x: C["CV_PEAK_UPPER"] <= x < C["CV_FLAT_UPPER"]
CV_VALL_CUTOFF = lambda x: C["CV_FLAT_UPPER"] <= x <= 1.00

def addSSHeader(pdbid, HELIX, SHEET):
    # Add HELIX/SHEET remarks to PDB header
    PDB_LINES = open("{}.pdb".format(pdbid)).readlines()
    PDB_OUT = open("{}.pdb".format(pdbid),'w')
    
    for j in range(len(HELIX)):
        PDB_OUT.write(HELIX[j]+"\n")
    
    for j in range(len(SHEET)):
        PDB_OUT.write(SHEET[j]+"\n")
    
    for j in range(len(PDB_LINES)):
        PDB_OUT.write(PDB_LINES[j])
    
    PDB_OUT.close()

def emptyDict(fields1, fields2=None, total1=False, total2=False, ensemble=False, N=1):
    data = {}
    if(fields2):
        for f1 in fields1:
            data[f1] = {}
            for f2 in fields2:
                if(ensemble):
                    data[f1][f2] = [0]*N
                else:
                    data[f1][f2] = 0
            if(total2):
                if(ensemble):
                    data[f1]['total'] = [0]*N
                else:
                    data[f1]['total'] = 0
    else:
        for f1 in fields1:
            if(ensemble):
                data[f1] = [0]*N
            else:
                data[f1] = 0
    if(total1):
        if(ensemble):
            data['total'] = [0]*N
        else:
            data['total'] = 0
    return data

def sumDicts(acc, val, sumTotalKeys=None, ensemble=False):
    for key in acc:
        if(key not in val):
            continue
        if(isinstance(acc[key], dict)):
            sumDicts(acc[key], val[key], sumTotalKeys=sumTotalKeys, ensemble=ensemble)
        else:
            if(ensemble):
                # Check if val is a list
                if(isinstance(val[key], list)):
                    for i in range(len(acc[key])):
                        acc[key][i] += val[key][i]
                    if(sumTotalKeys and (key in sumTotalKeys) and ('total' in acc)):
                        for i in range(len(acc['total'])):
                            acc['total'][i] += val[key][i]
                
                # Assume val is a float if not a list
                elif(isinstance(acc[key][index], list)):
                    for i in range(len(acc[key][index])):
                        acc[key][i] += val[key]
                    if(sumTotalKeys and (key in sumTotalKeys) and ('total' in acc)):
                        for i in range(len(acc['total'])):
                            acc['total'][i] += val[key]
                else:
                    acc[key][index] += val[key]
                    if(sumTotalKeys and (key in sumTotalKeys) and ('total' in acc)):
                        acc['total'][index] += val[key]
            else:
                acc[key] += val[key]
                if(sumTotalKeys and (key in sumTotalKeys) and ('total' in acc)):
                    acc['total'] += val[key]

def _getResiduePosition(residue):
    if("CA" in residue):
        return residue["CA"].get_coord()
    elif("N" in residue and "C" in residue):
        return (residue["N"].get_coord() + residue["C"].get_coord())/2
    elif("N" in residue and "O" in residue):
        return (residue["N"].get_coord() + residue["O"].get_coord())/2
    elif("C" in residue):
        return residue["C"].get_coord()
    else:
        raise Exception("Coordinate could not be defined for residue {}!".format(getID(residue=residue)))

def getSSECoordinate(model, sse):
    coord = np.zeros(3)
    for res_id in sse["residue_ids"]:
        ch, num, ins = res_id.split(".")
        rid = (' ', int(num), ins)
        residue = model[ch][rid]
        coord += _getResiduePosition(residue)
    return coord/len(sse["residue_ids"])

def getHelicoidalCoordinates(haxis, fixed_P, S, Gx, Gy, Gz, coord):
    m, n = haxis.shape
    mindist = (haxis[1,0]-coord[0])**2 + (haxis[1,1]-coord[1])**2 + (haxis[1,2]-coord[2])**2
    i_min = 0
    for i in range(2,m):
        d2 = (haxis[i,0]-coord[0])**2 + (haxis[i,1]-coord[1])**2 + (haxis[i,2]-coord[2])**2
        if(d2 < mindist):
            i_min = i
            mindist = d2
    
    # Get coordinates at minimum point
    sp = 3
    sm = 2
    if(i_min-2 < 0):
        sm = i_min
    if(i_min+3 > len(Gx)):
        sp =len(Gx)-i_min
    T = np.array([
        Gx[i_min-sm:i_min+sp].mean(),
        Gy[i_min-sm:i_min+sp].mean(),
        Gz[i_min-sm:i_min+sp].mean(),
    ])
    T /= np.linalg.norm(T) # tangent to helical axis at the point i_min
    
    # Check that angle between coord and T is within range (ideally pi/2)
    R = coord - haxis[i_min]
    h_angle = np.arccos(np.dot(R, T)/np.linalg.norm(R))
    if(h_angle < C["HELICOIDAL_ANGLE_CUTOFF_LOWER"] or h_angle > C["HELICOIDAL_ANGLE_CUTOFF_UPPER"]):
        # helicoidal coordinates undefined 
        return None
    
    # Compute helicoidal coordinates
    Xp = fixed_P - haxis[i_min]
    Xp -= np.dot(Xp,T)*T
    Xp /= np.linalg.norm(Xp)
    
    Yp = np.cross(T,Xp)
    phi_min = ma.atan2(np.dot(R,Yp), np.dot(R,Xp))
    
    return (round(phi_min,3), round(np.linalg.norm(R),3), round(S[i_min],3))

def getCM(model, rids):
    CM = np.array([0,0,0],dtype=float)
    M = 0
    for res_id in rids:
        cid, num, ins = res_id.split('.')
        rid = (' ', int(num), ins)
        if(rid in model[cid]):
            residue = model[cid][rid]
            for atom in residue:
                M += atom.mass
                CM += atom.mass*atom.get_coord()
    CM /= M
    return CM

def getPrincipalAxis(model, rids, CM):
    I = np.zeros((3,3), dtype=float)
    dmax = 0
    for rid in rids:
        cid, num, ins = rid.split('.')
        rid = (' ', int(num), ins)
        residue = model[cid][rid]
        for atom in residue:
            coord = atom.get_coord() - CM
            mass= atom.mass
            I[0,0] += mass*(coord[1]**2 + coord[2]**2)
            I[0,1] -= mass*coord[0]*coord[1]
            I[0,2] -= mass*coord[0]*coord[2]
            I[1,1] += mass*(coord[0]**2 + coord[2]**2)
            I[1,2] -= mass*coord[1]*coord[2]
            I[2,2] += mass*(coord[0]**2 + coord[1]**2)
            dmax = max(dmax, np.linalg.norm(coord))
    I[1,0] = I[0,1]
    I[2,0] = I[0,2]
    I[2,1] = I[1,2]
    w, v = np.linalg.eig(I)
    return w, v, dmax

def getFixedPoint(model, res_ids, haxis):
    CM = getCM(model, res_ids) # fixed point in space, determined by protein CM
    m, n = haxis.shape
    
    # Check if CM too close to DNA axis
    mindist = (haxis[0,0]-CM[0])**2 + (haxis[0,1]-CM[1])**2 + (haxis[0,2]-CM[2])**2
    i_min = 0
    for i in range(1,m):
        d2 = (haxis[i,0]-CM[0])**2 + (haxis[i,1]-CM[1])**2 + (haxis[i,2]-CM[2])**2
        if(d2 < mindist):
            mindist = d2
            i_min = i
    DCM = ma.sqrt(mindist)
    
    if(DCM < 1.0):
        # use protein principal axis to define fixed point
        w, v, _ = getPrincipalAxis(model, res_ids, CM)
        
        # get eigenvector corresponding to largest eigenvalue
        jmin = w.argmin()
        vmin = v[:,jmin]
        
        u = haxis[i_min] - CM
        u /= np.linalg.norm(u)
        dp = np.dot(u,vmin)
        if(dp >= 0):
            return CM + 20*vmin
        else:
            return CM - 20*vmin
    else:
        # otherwise use the center of mass
        return CM

def generateHaxis(haxis):
    # Generate haxis coordinates
    xcoef = np.array(haxis["x_coef"])
    ycoef = np.array(haxis["y_coef"])
    zcoef = np.array(haxis["z_coef"])
    
    m = int(haxis["axis_length"]*5) # 5 points per angstrom
    p = np.array(list(range(xcoef.size)))
    
    # Generate axis coordinates
    coords = np.zeros((m, 3))
    t = np.linspace(-haxis["axis_length"]/2, haxis["axis_length"]/2, num=m)
    for i in range(0, m):
        coords[i,0] = np.dot(xcoef, np.power(t[i], p))
        coords[i,1] = np.dot(ycoef, np.power(t[i], p))
        coords[i,2] = np.dot(zcoef, np.power(t[i], p))
    
    # Compute distance along axis
    S = np.zeros(m, dtype=float)
    for i in range(1,m):
        S[i] = S[i-1] + np.linalg.norm(coords[i]-coords[i-1])
    
    # Compute gradients of axis
    dS = np.diff(S)
    Gx = np.gradient(coords[1:,0], dS.mean())
    Gy = np.gradient(coords[1:,1], dS.mean())
    Gz = np.gradient(coords[1:,2], dS.mean())
    
    return {"coords": coords , "S": S, "Gx": Gx, "Gy": Gy, "Gz": Gz, "length": m}

def combineData(array, key, fields, dest_dict, delete=None):
    for item in array:
        if(delete):
            for f in delete:
                if(f in item):
                    del item[f]
        item_id = item[key]
        if(item_id not in dest_dict):
            dest_dict[item_id] = copy.deepcopy(item)
            for f in fields:
                dest_dict[item_id][f] = []
        
        # Loop over mutable fields
        for f in fields:
            dest_dict[item_id][f].append(item[f])

def getCVRatios(mi, res_ids, residues, key):
    ratios = {
        "peak_ratio": 0.0,
        "flat_ratio": 0.0,
        "valley_ratio": 0.0
    }
    srf_sesa = 0.0
    for rid in res_ids:
        if(not residues[rid]["surface"][mi]):
            continue
        cv = residues[rid][key][mi]
        if(cv is None):
            continue
        sesa = residues[rid]['sesa'][mi]["total"]
        srf_sesa += sesa
        if(CV_PEAK_CUTOFF(cv)):
            ratios["peak_ratio"] += sesa
        elif(CV_FLAT_CUTOFF(cv)):
            ratios["flat_ratio"] += sesa
        elif(CV_VALL_CUTOFF(cv)):
            ratios["valley_ratio"] += sesa
    ratios["peak_ratio"] /= srf_sesa
    ratios["flat_ratio"] /= srf_sesa
    ratios["valley_ratio"] /= srf_sesa
    
    return ratios

def getPropensities(mi, res_ids, residueBASA, residueFASA):
    propensity = {
        'ALA':None, 'CYS':None, 'ASP':None, 'GLU':None,
        'PHE':None, 'GLY':None, 'HIS':None, 'ILE':None,
        'LYS':None, 'LEU':None, 'MET':None, 'ASN':None,
        'PRO':None, 'GLN':None, 'ARG':None, 'SER':None,
        'THR':None, 'VAL':None, 'TRP':None, 'TYR':None
    }
    
    INT_FASA = {
        'ALA':0.0, 'CYS':0.0, 'ASP':0.0, 'GLU':0.0,
        'PHE':0.0, 'GLY':0.0, 'HIS':0.0, 'ILE':0.0,
        'LYS':0.0, 'LEU':0.0, 'MET':0.0, 'ASN':0.0,
        'PRO':0.0, 'GLN':0.0, 'ARG':0.0, 'SER':0.0,
        'THR':0.0, 'VAL':0.0, 'TRP':0.0, 'TYR':0.0
    }
    
    SRF_FASA = {
        'ALA':0.0, 'CYS':0.0, 'ASP':0.0, 'GLU':0.0,
        'PHE':0.0, 'GLY':0.0, 'HIS':0.0, 'ILE':0.0,
        'LYS':0.0, 'LEU':0.0, 'MET':0.0, 'ASN':0.0,
        'PRO':0.0, 'GLN':0.0, 'ARG':0.0, 'SER':0.0,
        'THR':0.0, 'VAL':0.0, 'TRP':0.0, 'TYR':0.0
    }
    
    int_fasa = 0.0
    srf_fasa = 0.0
    for rid in res_ids:
        resn = residueFASA[rid]['name']
        if(resn not in INT_FASA):
            continue
        if(not residueFASA[rid]["surface"][mi]):
            continue
        sc_fasa = residueFASA[rid]['fasa'][mi]['sc']
        if(rid in residueBASA):
            sc_basa = residueBASA[rid]['basa_sum']['sc']
        else:
            sc_basa = 0.0
        
        if(sc_basa > 0.0):
            INT_FASA[resn] += sc_fasa
            int_fasa += sc_fasa
        else:
            SRF_FASA[resn] += sc_fasa
            srf_fasa += sc_fasa
    if(srf_fasa == 0 or int_fasa == 0):
        # all or no residues interact with DNA - propensities are not defined
        return propensity 
    
    for key in propensity:
        if(INT_FASA[key] + SRF_FASA[key] == 0.0):
            # this residue doesn't appear in the interface or surface
            continue
        else:
            propensity[key] = round(np.log((1+INT_FASA[key]/int_fasa)/(1+SRF_FASA[key]/srf_fasa)),4)
    
    return propensity

def _removeItem(mi, DATA, keys):
    if(len(keys) == 0):
        DATA.pop(mi)
    elif(isinstance(DATA, list)):
        key = keys.pop(0)
        for D in DATA:
            _removeItem(mi, D[key], keys)
    else:
        _removeItem(mi, DATA[keys.pop(0)], keys)

def deleteModelData(mi, DATA, fields):
    for field in fields:
        keys = field.split('.')
        _removeItem(mi, DATA[keys.pop(0)], keys)

def comp(pdbid, N, assembly, PRO_DATA, DSSP, DNA_DATA, INT_DATA, REGEXES, NUCLEOTIDES, COMPONENTS, META_DATA,
    mmcif_dict=None,
    ADD_MMCIF=False,
    REMOVE_HEADER=True
):
    #----------------------------------Template Dictionaries---------------------------------------#
    JSON = {
        "structure_id": pdbid,
        "num_models": N,
        "protein": {
            "chains": {},
            "residues": {},
            "num_chains": None,
            "num_residues": None,
            "models": []
        },
        "dna": {
            "nucleotides": {},
            "chains": None,
            "num_chains": None,
            "num_nucleotides": None,
            "models": []
        },
        "interfaces": {
            "models": []
        },
        "meta_data": {}
    }
    
    interface_template = {
        "dna_entity_id": None,
        "pro_entity_id": None,
        "protein_chains": set(),
        "nucleotide-residue_interactions": {},
        "nucleotide-sse_interactions": {},
        "interface_features": {},
        "nucleotide_data": {},
        "residue_data": {},
        "sse_data": {}
    }
    
    interface_features_template = {
        "basa": emptyDict(nucMtyLabel, fields2=resMtyLabel, total1=True),
        "hbond_sum": emptyDict(nucMtyLabel, fields2=resMtyLabel, total1=True),
        "vdw_sum": emptyDict(nucMtyLabel, fields2=resMtyLabel, total1=True),
        "secondary_structure_composition": None,
        "residue_propensities": {},
        "protein_surface_geometry": {
            "cv_fine": None,
            "cv_coarse": None
        },
        "mean_hydrophobicity_score": 0.0,
        "psuedo-stack_interaction_ratio": None,
        "psuedo-pair_interaction_ratio": None,
        "protein_chain_id": None,
        "interaction_count": 0,
        "weak_interaction_count": 0,
        "stack_count": 0,
        "pair_count": 0,
        "weight_sum": 0.0,
        "residue_ids": set(),
        "segment_ids": set(),
        "interaction_moiety_summary": emptyDict(nucMtyLabel, fields2=resSST)
    }
    interface_features_template["basa"]["secondary_structure"] =  emptyDict(nucMtyLabel, fields2=resSST)
    interface_features_template["hbond_sum"]["secondary_structure"] =  emptyDict(nucMtyLabel, fields2=resSST)
    interface_features_template["vdw_sum"]["secondary_structure"] =  emptyDict(nucMtyLabel, fields2=resSST)
    
    interaction_template = {
        "basa": None,
        "nuc_basa": None,
        "nuc_chain": None,
        "nuc_id": None,
        "nuc_name": None,
        "nuc_number": None,
        "res_basa": None,
        "res_chain": None,
        "res_id": None,
        "res_name": None,
        "res_number": None,
        "res_secondary_structure": None,
        "nuc_secondary_structure": None,
        "geometry": 'none',
        "hbonds": [],
        "vdw_interactions": [],
        "hbond_sum": None,
        "vdw_sum": None,
        "min_distance": None,
        "mean_nn_distance": None,
        "weak_interaction": False
    }
    
    sse_interaction_template = {
        "basa": None,
        "nuc_basa": None,
        "nuc_chain": None,
        "nuc_id": None,
        "nuc_name": None,
        "nuc_number": None,
        "sse_basa": None,
        "sse_chain": None,
        "sse_id": None,
        "hbond_sum": None,
        "vdw_sum": None,
    }
    
    residue_sum_template = {
        "res_id": None,
        "basa_sum": emptyDict(resMtyLabel, total1=True),
        "hbond_sum": emptyDict(nucMtyLabel, fields2=resMtyLabel, total1=True),
        "vdw_interaction_sum": emptyDict(nucMtyLabel, fields2=resMtyLabel, total1=True),
        "nucleotide_interaction_count": 0,
        "interacts_with": set(),
        "interacts_by": set()
    }
    
    nucleotide_sum_template = {
        "nuc_id": None,
        "basa_sum": None,
        "hbond_sum": None,
        "vdw_interaction_sum": None,
        "residue_interaction_count": 0,
        "interacts_with": set(),
        "interacts_by": set()
    }
    #----------------------------------Template Dictionaries---------------------------------------#
    
    ### Consolidate immutable data from protein data and append mutables ###
    chain_mutables = ["secondary_structure", "continuous", "interacts_with_dna"]
    residue_mutables = ["secondary_structure", "sap_score", "cv_fine", "cv_coarse", "sesa", "surface"]
    pro_delete = ["chains", "residues"]
    for model in PRO_DATA:
        # Add chain data
        combineData(model["chains"], "id", chain_mutables, JSON["protein"]["chains"])
        
        # Add residue data
        combineData(model["residues"], "id", residue_mutables, JSON["protein"]["residues"])
        
        # Delete certain fields from model
        for key in pro_delete:
            del model[key]
        JSON["protein"]["models"].append(model)
    
    ### Add nmer-residue lookup and residue-sse lookup info ###
    res_sse_lookup = [] # a list of dicts, one for each model, that maps resid to parent SSE
    pro_entity_lookup = [] # a list of dicts, one for each model, that maps resid to parent entity
    res_segment_lookup = [] # a list of dicts, one for each model, that maps resid to parent segment
    pro_entity_map = [] # a list of dicts, one for each model, that maps entity ID to entity
    for model in JSON["protein"]["models"]:
        # map resid to parent sse
        sse_lookup = {}
        res_lookup = {}
        eid_lookup = {}
        seg_lookup = {}
        for sse in model["secondary_structure_elements"]:
            for res_id in sse["residue_ids"]:
                sse_lookup[res_id] = sse
        for entity in model["entities"]:
            eid = entity["id"]
            for res_id in entity["residue_ids"]:
                res_lookup[res_id] = eid
            eid_lookup[eid] = entity
        for seg in model["segments"]:
            for res_id in seg["residue_ids"]:
                seg_lookup[res_id] = seg["id"]
        res_sse_lookup.append(sse_lookup)
        pro_entity_lookup.append(res_lookup)
        pro_entity_map.append(eid_lookup)
        res_segment_lookup.append(seg_lookup)
        
    JSON["protein"]["num_chains"] = len(JSON["protein"]["chains"])
    JSON["protein"]["num_residues"] = len(JSON["protein"]["residues"])
    
    ### Consolidate immutable data from dna data and append mutables ###
    JSON["dna"]["chains"] = DNA_DATA[0]["chains"]
    JSON["dna"]["num_nucleotides"] = DNA_DATA[0]["num_nucleotides"]
    JSON["dna"]["num_chains"] = DNA_DATA[0]["num_chains"]
    nucleotide_mutables = ["glycosidic_conformation", "origin", "secondary_structure", "graph_coordinates"]
    dna_delete = ["num_nucleotides", "num_chains", "chains", "nucleotides"]
    for model in DNA_DATA:
        # Add nucleotide data
        combineData(model["nucleotides"], "id", nucleotide_mutables, JSON["dna"]["nucleotides"], delete=["groove_atoms"])
        
        # Delete certain fields from model
        for key in dna_delete:
            del model[key]
        JSON["dna"]["models"].append(model)
    
    ### Add entity-nucleotide lookup info ###
    dna_entity_lookup = [] # a list of dicts, one for each model which maps nucleotide id to parent entity id
    dna_entity_map = [] # a list of dicts, one for each model which maps entity id to entity
    for model in JSON["dna"]["models"]:
        nuc_lookup = {}
        eid_lookup = {}
        for entity in model["entities"]:
            eid = entity["id"]
            eid_lookup[eid] = entity
            for n in entity["nucleotides"]:
                nuc_lookup[n] = eid
        dna_entity_lookup.append(nuc_lookup)
        dna_entity_map.append(eid_lookup)
    
    ### Make nucleotide map ###
    #nucleotideMap = {}
    #for nuc in JSON["dna"]["nucleotides"].values():
    #    nucleotideMap[nuc["id"]] = nuc
    
    ### Add nucleotide and residue FASA values from INT_DATA ###
    for mi in range(N):
        model = INT_DATA[mi]
        
        for res in model["basa"]["residues"]:
            rid = res["id"]
            fasa = res["fasa"]
            if(JSON["protein"]["residues"][rid]["fasa"] is None):
                JSON["protein"]["residues"][rid]["fasa"] = []
            JSON["protein"]["residues"][rid]["fasa"].append(fasa)
        
        for nuc in model["basa"]["nucleotides"]:
            nid = nuc["id"]
            fasa = nuc["fasa"]
            if(JSON["dna"]["nucleotides"][nid]["fasa"] is None):
                JSON["dna"]["nucleotides"][nid]["fasa"] = []
            JSON["dna"]["nucleotides"][nid]["fasa"].append(fasa)
    
    DELETE_MODELS = []
    ### Add nucleotide-residue interactions ###
    for mi in range(N):
        model = INT_DATA[mi]
        interfaces = {} # keyed by entity_id
        
        ### Add interactions, starting from interaction list ###
        for nr in model["nucleotide-residue_interactions"]:
            # Check if in valid DNA entity
            if(not (nr["nuc_id"] in dna_entity_lookup[mi])):
                continue
            if(not (nr["res_id"] in pro_entity_lookup[mi])):
                continue
            
            # Check if this is a surface residue
            if(not JSON["protein"]["residues"][nr["res_id"]]["surface"][mi]):
                continue
            
            res_chain = nr["res_id"][0]
            nr_id = getHash(nr["nuc_id"], nr["res_id"])
            int_id = getHash(pro_entity_lookup[mi][nr["res_id"]], dna_entity_lookup[mi][nr["nuc_id"]])
            if(int_id not in interfaces):
                # Add a new interface
                interfaces[int_id] = copy.deepcopy(interface_template)
                interfaces[int_id]["dna_entity_id"] = dna_entity_lookup[mi][nr["nuc_id"]]
                interfaces[int_id]["pro_entity_id"] = pro_entity_lookup[mi][nr["res_id"]]
            interfaces[int_id]["protein_chains"].add(res_chain)
            
            # Add interaction to interface
            interfaces[int_id]["nucleotide-residue_interactions"][nr_id] = copy.deepcopy(interaction_template)
            interfaces[int_id]["nucleotide-residue_interactions"][nr_id].update(nr)
            
        # Add BASA information
        for nr in model["basa"]["interactions"]:
            if(not (nr["nuc_id"] in dna_entity_lookup[mi])):
                continue
            if(not (nr["res_id"] in pro_entity_lookup[mi])):
                continue
            
            nr_id = getHash(nr["nuc_id"], nr["res_id"])
            int_id = getHash(pro_entity_lookup[mi][nr["res_id"]], dna_entity_lookup[mi][nr["nuc_id"]])
            if(int_id in interfaces):
                if(nr_id in interfaces[int_id]["nucleotide-residue_interactions"]):
                    interfaces[int_id]["nucleotide-residue_interactions"][nr_id].update(nr)
        
        # Add hydrogen bonds
        for hb in model["hbond"]:
            # check if in valid DNA entity
            if(not (hb["nuc_id"] in dna_entity_lookup[mi])):
                continue
            if(not (hb["res_id"] in pro_entity_lookup[mi])):
                continue
            
            nr_id = getHash(hb["nuc_id"], hb["res_id"])
            int_id = getHash(pro_entity_lookup[mi][hb["res_id"]], dna_entity_lookup[mi][hb["nuc_id"]])
            # Add hbond to interaction
            if(int_id in interfaces):
                if(nr_id in interfaces[int_id]["nucleotide-residue_interactions"]):
                    interfaces[int_id]["nucleotide-residue_interactions"][nr_id]["hbonds"].append(
                        {
                            "res_atom": hb["res_atom"],
                            "distance": hb["distance"],
                            "nuc_atom": hb["nuc_atom"],
                            "nuc_moiety": hb["nuc_moiety"],
                            "res_moiety": hb["res_moiety"]
                        }
                    )
        
        # Add vdw interactions
        for vdw in model["vdw"]:
            # check if in valid DNA entity
            if(not (vdw["nuc_id"] in dna_entity_lookup[mi])):
                continue
            if(not (vdw["res_id"] in pro_entity_lookup[mi])):
                continue
            
            nr_id = getHash(vdw["nuc_id"], vdw["res_id"])
            int_id = getHash(pro_entity_lookup[mi][vdw["res_id"]], dna_entity_lookup[mi][vdw["nuc_id"]])
            # Add vdw to interaction
            if(int_id in interfaces):
                if(nr_id in interfaces[int_id]["nucleotide-residue_interactions"]):
                    interfaces[int_id]["nucleotide-residue_interactions"][nr_id]["vdw_interactions"].append(
                        {
                            "res_atom": vdw["res_atom"],
                            "distance": vdw["distance"],
                            "nuc_atom": vdw["nuc_atom"],
                            "nuc_moiety": vdw["nuc_moiety"],
                            "res_moiety": vdw["res_moiety"]
                        }
                    )
        
        # Add geometries
        for geo in model["geometry"]:
            # check if in valid DNA entity
            if(not (geo["nuc_id"] in dna_entity_lookup[mi])):
                continue
            if(not (geo["res_id"] in pro_entity_lookup[mi])):
                continue
                
            nr_id = getHash(geo["nuc_id"], geo["res_id"])
            int_id = getHash(pro_entity_lookup[mi][geo["res_id"]], dna_entity_lookup[mi][geo["nuc_id"]])
            # Add geometry to interaction
            if(int_id in interfaces):
                if(nr_id in interfaces[int_id]["nucleotide-residue_interactions"]):
                    interfaces[int_id]["nucleotide-residue_interactions"][nr_id]["geometry"] = geo["geometry"]
        
        ### Compile interface statistics ###
        for int_id in interfaces:
            # Add interface features dicts
            for chain in interfaces[int_id]["protein_chains"]:
                interfaces[int_id]["interface_features"][chain] = copy.deepcopy(interface_features_template)
                interfaces[int_id]["interface_features"][chain]["protein_chain_id"] = chain
            
            # Compile interaction statistics over each interface
            sse_nuc_count = {}
            for nr_id in interfaces[int_id]["nucleotide-residue_interactions"]:
                nr = interfaces[int_id]["nucleotide-residue_interactions"][nr_id]
                nuc_id = nr["nuc_id"]
                res_id = nr["res_id"]
                parent_sse_id = res_sse_lookup[mi][res_id]["id"]
                segment_id = res_segment_lookup[mi][res_id]
                res_chain = res_id[0]
                
                interfaces[int_id]["interface_features"][res_chain]["residue_ids"].add(res_id)
                interfaces[int_id]["interface_features"][res_chain]["segment_ids"].add(segment_id)
                nuc_ss = JSON["dna"]["nucleotides"][nuc_id]["secondary_structure"][mi]
                res_ss = JSON["protein"]["residues"][res_id]["secondary_structure"][mi]
                if(nuc_ss == "helical"):
                    nucMty = nucDSMtyLabel
                else:
                    nucMty = nucSSMtyLabel
                nr["hbond_sum"] = emptyDict(nucMty, fields2=resMtyLabel, total1=True)
                nr["vdw_sum"] = emptyDict(nucMty, fields2=resMtyLabel, total1=True)
                
                # Add Secondary Structure
                nr["res_secondary_structure"] = res_ss
                nr["nuc_secondary_structure"] = nuc_ss
                
                # Add empty basa basa dictionaries if needed
                if(nr["basa"] is None):
                    nr["basa"] = emptyDict(nucMty, fields2=resMtyLabel, total1=True)
                    nr["res_basa"] = emptyDict(resMtyLabel, total1=True)
                    nr["nuc_basa"] = emptyDict(nucMty, total1=True)
                    
                # Sum hydrogen bonds
                for hb in nr["hbonds"]:
                    nm = hb["nuc_moiety"]
                    rm = hb["res_moiety"]
                    if(rm and nm):
                        nr["hbond_sum"][nm][rm] += 1
                        nr["hbond_sum"]["total"] += 1
                
                # Sum vdw interactions
                for vw in nr["vdw_interactions"]:
                    nm = vw["nuc_moiety"]
                    rm = vw["res_moiety"]
                    if(rm and nm):
                        nr["vdw_sum"][nm][rm] += 1
                        nr["vdw_sum"]["total"] += 1
                
                # Classify interaction moieties
                nname = nr["nuc_name"]
                rname = nr["res_name"]
                if(not re.search(REGEXES["PROTEIN"]["STANDARD_RESIDUES"], rname)):
                    rname = COMPONENTS[rname]['_chem_comp.mon_nstd_parent_comp_id']
                if(not re.search(REGEXES["DNA"]["STANDARD_NUCLEOTIDES"], nname)):
                    nname = COMPONENTS[nname]['_chem_comp.mon_nstd_parent_comp_id']
                
                if(INTERACTION_STATS[nname+rname]["count"] >= 100):
                    # Classify interaction
                    int_mty, nuc_int_mty, res_int_mty = getInteractionMoiety(nr, INTERACTION_STATS[nname+rname], ["basa", "hbond_sum", "vdw_sum"])
                else:
                    # Try to classify interactions with fallback stats.
                    int_mty, nuc_int_mty, res_int_mty = getInteractionMoiety(nr, INTERACTION_STATS["fallback"], ["basa", "hbond_sum", "vdw_sum"])
                if(len(int_mty) == 0):
                    # A weak interaction - flag it, define interaction moiety based on max basa
                    nr["weak_interaction"] = True
                    max_basa = 0.0
                    max_nm = None
                    max_rm = None
                    for nm in nr['basa']:
                        if(nm == "total"):
                            continue
                        for rm in nr['basa'][nm]:
                            if(nr['basa'][nm][rm] > max_basa):
                                max_basa = nr['basa'][nm][rm]
                                max_nm = nm
                                max_rm = rm
                    int_mty = ["{}.{}".format(max_nm, max_rm)]
                    nuc_int_mty = [max_nm]
                    res_int_mty = [max_rm]
                nr["moiety_interactions"] = int_mty
                nr["nucleotide_interaction_moieties"] = nuc_int_mty
                nr["residue_interaction_moieties"] = res_int_mty
                
                # Add to interaction count
                interfaces[int_id]["interface_features"][res_chain]["interaction_count"] += 1
                if(nr["weak_interaction"]):
                    interfaces[int_id]["interface_features"][res_chain]["weak_interaction_count"] += 1
                
                # Add to moiety interaction summary
                if(not nr["weak_interaction"]):
                    for m in nr["nucleotide_interaction_moieties"]:
                        interfaces[int_id]["interface_features"][res_chain]["interaction_moiety_summary"][m][nr["res_secondary_structure"]] += 1
                
                # Sum basa, vdw, hbonds
                sumDicts(interfaces[int_id]["interface_features"][res_chain]["basa"], nr["basa"], sumTotalKeys=resMtyLabel)
                sumDicts(interfaces[int_id]["interface_features"][res_chain]["hbond_sum"], nr["hbond_sum"])
                sumDicts(interfaces[int_id]["interface_features"][res_chain]["vdw_sum"], nr["vdw_sum"])
                
                # Add secondary structure sums
                for mn in nucMty:
                    for mr in resMtyLabel:
                        interfaces[int_id]["interface_features"][res_chain]["basa"]["secondary_structure"][mn][res_ss] += nr["basa"][mn][mr]
                        interfaces[int_id]["interface_features"][res_chain]["hbond_sum"]["secondary_structure"][mn][res_ss] += nr["hbond_sum"][mn][mr]
                        interfaces[int_id]["interface_features"][res_chain]["vdw_sum"]["secondary_structure"][mn][res_ss] += nr["vdw_sum"][mn][mr]
                
                # sum nucleotide properties
                if(nuc_id not in interfaces[int_id]["nucleotide_data"]):
                    # Make new entry
                    interfaces[int_id]["nucleotide_data"][nuc_id] = copy.deepcopy(nucleotide_sum_template)
                    interfaces[int_id]["nucleotide_data"][nuc_id]["nuc_id"] = nuc_id
                    interfaces[int_id]["nucleotide_data"][nuc_id]["basa_sum"] =  emptyDict(nucMty, total1=True)
                    interfaces[int_id]["nucleotide_data"][nuc_id]["basa_sum"]["secondary_structure"] = emptyDict(nucMty, fields2=resSST)
                    interfaces[int_id]["nucleotide_data"][nuc_id]["hbond_sum"] =  emptyDict(nucMty, fields2=resMtyLabel, total1=True)
                    interfaces[int_id]["nucleotide_data"][nuc_id]["hbond_sum"]["secondary_structure"] = emptyDict(nucMty, fields2=resSST)
                    interfaces[int_id]["nucleotide_data"][nuc_id]["vdw_interaction_sum"] =  emptyDict(nucMty, fields2=resMtyLabel, total1=True)
                    interfaces[int_id]["nucleotide_data"][nuc_id]["vdw_interaction_sum"]["secondary_structure"] = emptyDict(nucMty, fields2=resSST)
                
                # Update various accumulator dictionaries
                interfaces[int_id]["nucleotide_data"][nuc_id]["residue_interaction_count"] += 1
                sumDicts(interfaces[int_id]["nucleotide_data"][nuc_id]["basa_sum"], nr["nuc_basa"])
                sumDicts(interfaces[int_id]["nucleotide_data"][nuc_id]["hbond_sum"], nr["hbond_sum"], sumTotalKeys=resMtyLabel)
                sumDicts(interfaces[int_id]["nucleotide_data"][nuc_id]["vdw_interaction_sum"], nr["vdw_sum"], sumTotalKeys=resMtyLabel)
                interfaces[int_id]["nucleotide_data"][nuc_id]["interacts_with"].update(res_int_mty)
                interfaces[int_id]["nucleotide_data"][nuc_id]["interacts_by"].update(nuc_int_mty)
                
                # Add secondary structure sums
                for mn in nucMty:
                    interfaces[int_id]["nucleotide_data"][nuc_id]["basa_sum"]["secondary_structure"][mn][res_ss] += nr["nuc_basa"][mn]
                    for mr in resMtyLabel:
                        interfaces[int_id]["nucleotide_data"][nuc_id]["hbond_sum"]["secondary_structure"][mn][res_ss] += nr["hbond_sum"][mn][mr]
                        interfaces[int_id]["nucleotide_data"][nuc_id]["vdw_interaction_sum"]["secondary_structure"][mn][res_ss] += nr["hbond_sum"][mn][mr]
                
                # Sum residue properties
                if(res_id not in interfaces[int_id]["residue_data"]):
                    # Make new entry
                    interfaces[int_id]["residue_data"][res_id] = copy.deepcopy(residue_sum_template)
                    interfaces[int_id]["residue_data"][res_id]["res_id"] = res_id
                    interfaces[int_id]["residue_data"][res_id]["interacting_nucleotides"] = set()
                
                # Update various accumulator dictionaries
                interfaces[int_id]["residue_data"][res_id]["nucleotide_interaction_count"] += 1
                sumDicts(interfaces[int_id]["residue_data"][res_id]["basa_sum"], nr["res_basa"])
                sumDicts(interfaces[int_id]["residue_data"][res_id]["hbond_sum"], nr["hbond_sum"], sumTotalKeys=resMtyLabel)
                sumDicts(interfaces[int_id]["residue_data"][res_id]["vdw_interaction_sum"], nr["vdw_sum"], sumTotalKeys=resMtyLabel)
                interfaces[int_id]["residue_data"][res_id]["interacts_with"].update(nuc_int_mty)
                interfaces[int_id]["residue_data"][res_id]["interacts_by"].update(res_int_mty)
                interfaces[int_id]["residue_data"][res_id]["interacting_nucleotides"].add(nuc_id)
                
                # Update geometry counts
                if(nr["geometry"] == "pseudo_stack"):
                    interfaces[int_id]["interface_features"][res_chain]["stack_count"] += 1.0
                elif(nr["geometry"] == "pseudo_pair"):
                    interfaces[int_id]["interface_features"][res_chain]["pair_count"] += 1.0
                
                # Check for existing SSE-nucleotide interactions
                ns_id = getHash(nuc_id, parent_sse_id)
                if(parent_sse_id not in sse_nuc_count):
                    sse_nuc_count[parent_sse_id] = 0
                if(ns_id not in interfaces[int_id]["nucleotide-sse_interactions"]):
                    # Make new entry
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id] = copy.deepcopy(sse_interaction_template)
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["sse_id"] = parent_sse_id
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["sse_chain"] = nr["res_chain"]
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["nuc_chain"] = nr["nuc_chain"]
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["nuc_name"] = nr["nuc_name"]
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["nuc_number"] = nr["nuc_number"]
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["nuc_id"] = nr["nuc_id"]
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["basa"] = emptyDict(nucMty, fields2=resMtyLabel, total1=True)
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["sse_basa"] = emptyDict(resMtyLabel, total1=True)
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["nuc_basa"] = emptyDict(nucMty, total1=True)
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["hbond_sum"] = emptyDict(nucMty, fields2=resMtyLabel, total1=True)
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["vdw_sum"] = emptyDict(nucMty, fields2=resMtyLabel, total1=True)
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["nucleotide_interaction_moieties"] = set()
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["residue_interaction_moieties"] = set()
                    interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["moiety_interactions"] = set()
                sumDicts(interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["sse_basa"], nr["res_basa"])
                sumDicts(interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["nuc_basa"], nr["nuc_basa"])
                sumDicts(interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["basa"], nr["basa"])
                sumDicts(interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["hbond_sum"], nr["hbond_sum"])
                sumDicts(interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["vdw_sum"], nr["vdw_sum"])
                interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["nucleotide_interaction_moieties"].update(nuc_int_mty)
                interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["residue_interaction_moieties"].update(res_int_mty)
                interfaces[int_id]["nucleotide-sse_interactions"][ns_id]["moiety_interactions"].update(int_mty)
                sse_nuc_count[parent_sse_id] += 1
                
                # Add interactions to parent strand/chain
                parent_dna_entity = dna_entity_map[mi][dna_entity_lookup[mi][nuc_id]]
                for s in parent_dna_entity["strands"]:
                    if(nuc_id in s["ids"]):
                        s["interacts_with_protein"] = True
                        break
                JSON["protein"]["chains"][res_chain]["interacts_with_dna"][mi] = True
            
            # Compute residue propensities and surface geometries
            for chain in interfaces[int_id]["protein_chains"]:
                res_ids = JSON["protein"]["chains"][chain]["residue_ids"]
                interfaces[int_id]["interface_features"][chain]["residue_propensities"] = getPropensities(mi, res_ids, interfaces[int_id]["residue_data"], JSON["protein"]["residues"])
                interfaces[int_id]["interface_features"][chain]["protein_surface_geometry"]["cv_fine"] = getCVRatios(mi, res_ids, JSON["protein"]["residues"], "cv_fine")
                interfaces[int_id]["interface_features"][chain]["protein_surface_geometry"]["cv_coarse"] = getCVRatios(mi, res_ids, JSON["protein"]["residues"], "cv_coarse")
            
            # Compile SSE data from residue data
            for rid in interfaces[int_id]["residue_data"]:
                parent_sse_id = res_sse_lookup[mi][rid]["id"]
                if(parent_sse_id not in interfaces[int_id]["sse_data"]):
                    # Add new entry
                    interfaces[int_id]["sse_data"][parent_sse_id] = copy.deepcopy(residue_sum_template)
                    interfaces[int_id]["sse_data"][parent_sse_id]["sse_id"] = parent_sse_id
                    interfaces[int_id]["sse_data"][parent_sse_id]["nucleotide_interaction_count"] = sse_nuc_count[parent_sse_id]
                    interfaces[int_id]["sse_data"][parent_sse_id]["interacting_nucleotides"] = set()
                    
                    # remove unused fields
                    del interfaces[int_id]["sse_data"][parent_sse_id]["res_id"]
                
                # Sum over all accumulated values
                sumDicts(interfaces[int_id]["sse_data"][parent_sse_id]["basa_sum"], interfaces[int_id]["residue_data"][rid]["basa_sum"])
                sumDicts(interfaces[int_id]["sse_data"][parent_sse_id]["hbond_sum"], interfaces[int_id]["residue_data"][rid]["hbond_sum"])
                sumDicts(interfaces[int_id]["sse_data"][parent_sse_id]["vdw_interaction_sum"], interfaces[int_id]["residue_data"][rid]["vdw_interaction_sum"])
                interfaces[int_id]["sse_data"][parent_sse_id]["interacts_with"].update(interfaces[int_id]["residue_data"][rid]["interacts_with"])
                interfaces[int_id]["sse_data"][parent_sse_id]["interacts_by"].update(interfaces[int_id]["residue_data"][rid]["interacts_by"])
                interfaces[int_id]["sse_data"][parent_sse_id]["interacting_nucleotides"].update(interfaces[int_id]["residue_data"][rid]["interacting_nucleotides"])
                
                
            # Geometry ratios
            for chain in interfaces[int_id]["protein_chains"]:
                interfaces[int_id]["interface_features"][chain]["psuedo-stack_interaction_ratio"] = (
                    interfaces[int_id]["interface_features"][chain]["stack_count"]/
                    interfaces[int_id]["interface_features"][chain]["interaction_count"]
                )
                interfaces[int_id]["interface_features"][chain]["psuedo-pair_interaction_ratio"] = (
                    interfaces[int_id]["interface_features"][chain]["pair_count"]/
                    interfaces[int_id]["interface_features"][chain]["interaction_count"]
                )
        
        ### Add helicoidal coordinates to each interface and update interaction moieties ###
        for int_id in interfaces:
            eid = interfaces[int_id]["dna_entity_id"]
            entity = dna_entity_map[mi][eid]
            for i in range(len(entity["helical_segments"])):
                # Compute helicoidal coordinates for SSE, residues and nucleotides
                HELIX = entity["helical_segments"][i]
                HAXIS = generateHaxis(HELIX["helical_axis"])
                FIXED = getFixedPoint(assembly[mi] , res_ids, HAXIS["coords"])
                nucleotide_ids = set(HELIX["ids1"] + HELIX["ids2"])
                
                # Residue coordinates
                for res_id in interfaces[int_id]["residue_data"]:
                    if(len(interfaces[int_id]["residue_data"][res_id]["interacting_nucleotides"] & nucleotide_ids) == 0):
                        # this residue doesn't interact with helix nucleotides
                        continue 
                    ch, num, ins = res_id.split('.')
                    rid = (' ', int(num), ins)
                    res = assembly[mi][ch][rid]
                    coord = _getResiduePosition(res)
                    hc = getHelicoidalCoordinates(HAXIS["coords"], FIXED, HAXIS["S"], HAXIS["Gx"], HAXIS["Gy"], HAXIS["Gz"], coord)
                    if(hc is None):
                        # skip residues which are too far away from the helical axis ends - not well defined 
                        continue
                    interfaces[int_id]["residue_data"][res_id]["helicoidal_coordinates"] = {
                        "phi": hc[0],
                        "rho": hc[1],
                        "s": hc[2]
                    }
                
                # SSE coordinates
                for sse in JSON["protein"]["models"][mi]["secondary_structure_elements"]:
                    sse_id = sse["id"]
                    if(sse_id in interfaces[int_id]["sse_data"]):
                        if(len(interfaces[int_id]["sse_data"][sse_id]["interacting_nucleotides"] & nucleotide_ids) == 0):
                            # this sse doesn't interact with helix nucleotides
                            continue 
                        coord = getSSECoordinate(assembly[mi], sse)
                        hc = getHelicoidalCoordinates(HAXIS["coords"], FIXED, HAXIS["S"], HAXIS["Gx"], HAXIS["Gy"], HAXIS["Gz"], coord)
                        if(hc is None):
                            # skip sse which are too far away from the helical axis ends - not well defined 
                            continue
                        interfaces[int_id]["sse_data"][sse_id]["helicoidal_coordinates"] = {
                            "phi": hc[0],
                            "rho": hc[1],
                            "s": hc[2]
                        }
                
                # Get binding sites
                bindsite1 = []
                bindsite2 = []
                for i in range(HELIX["length"]):
                    id1 = HELIX["ids1"][i]
                    id2 = HELIX["ids2"][i]
                    if(id1 in interfaces[int_id]["nucleotide_data"]):
                            bindsite1.append(i)
                    if(id2 in interfaces[int_id]["nucleotide_data"]):
                            bindsite2.append(i)
                
                if(len(bindsite1) > 0):
                    bsseq1 = HELIX["sequence1"][max(bindsite1[0]-2,0):min(bindsite1[-1]+3,len(HELIX["sequence1"]))]
                else:
                    bsseq1 = ''
                if(len(bindsite2) > 0):
                    bsseq2 = HELIX["sequence2"][max(bindsite2[0]-2,0):min(bindsite2[-1]+3,len(HELIX["sequence2"]))]
                else:
                    bsseq2 = ''
                # Add to interface features
                interfaces[int_id]["binding_site1"] = bsseq1
                interfaces[int_id]["binding_site2"] = bsseq2
            
            # Nucleotide interaction classifications
            for nuc_id in interfaces[int_id]["nucleotide_data"]:
                interfaces[int_id]["nucleotide_data"][nuc_id]["interacts_by"] = list(interfaces[int_id]["nucleotide_data"][nuc_id]["interacts_by"])
                interfaces[int_id]["nucleotide_data"][nuc_id]["interacts_with"] = list(interfaces[int_id]["nucleotide_data"][nuc_id]["interacts_with"])
            
            # Residue interaction classifications
            for res_id in interfaces[int_id]["residue_data"]:
                res_chain = res_id[0]
                interfaces[int_id]["residue_data"][res_id]["interacts_by"] = list(interfaces[int_id]["residue_data"][res_id]["interacts_by"])
                interfaces[int_id]["residue_data"][res_id]["interacts_with"] = list(interfaces[int_id]["residue_data"][res_id]["interacts_with"])
                interfaces[int_id]["residue_data"][res_id]["interacting_nucleotides"] = list(interfaces[int_id]["residue_data"][res_id]["interacting_nucleotides"])
                if(JSON["protein"]["residues"][res_id]["sap_score"][mi] is not None):
                    interfaces[int_id]["interface_features"][res_chain]["mean_hydrophobicity_score"] += (JSON["protein"]["residues"][res_id]["fasa"][mi]["sc"]
                        *JSON["protein"]["residues"][res_id]["sap_score"][mi])
                    interfaces[int_id]["interface_features"][res_chain]["weight_sum"] += JSON["protein"]["residues"][res_id]["fasa"][mi]["sc"]
            
            # SSE interaction classifications
            for sse_id in interfaces[int_id]["sse_data"]:
                interfaces[int_id]["sse_data"][sse_id]["interacts_with"] = list(interfaces[int_id]["sse_data"][sse_id]["interacts_with"])
                interfaces[int_id]["sse_data"][sse_id]["interacts_by"] = list(interfaces[int_id]["sse_data"][sse_id]["interacts_by"])
                interfaces[int_id]["sse_data"][sse_id]["interacting_nucleotides"] = list(interfaces[int_id]["sse_data"][sse_id]["interacting_nucleotides"])
            
            # Classify secondary structure composition
            for chain in interfaces[int_id]["protein_chains"]:
                interfaces[int_id]["interface_features"][chain]["residue_ids"] = list(interfaces[int_id]["interface_features"][chain]["residue_ids"])
                interfaces[int_id]["interface_features"][chain]["segment_ids"] = list(interfaces[int_id]["interface_features"][chain]["segment_ids"])
                totFASA = 0.0
                Hsum = 0.0
                Ssum = 0.0
                Lsum = 0.0
                
                for resID in interfaces[int_id]["interface_features"][chain]["residue_ids"]:
                    fasa = JSON["protein"]["residues"][resID]["fasa"][mi]["sc"]
                    totFASA += fasa
                    if(JSON["protein"]["residues"][resID]["secondary_structure"][mi] == 'H'):
                        Hsum += fasa
                    elif(JSON["protein"]["residues"][resID]["secondary_structure"][mi] == 'S'):
                        Ssum += fasa
                    elif(JSON["protein"]["residues"][resID]["secondary_structure"][mi] == 'L'):
                        Lsum += fasa
                Hsum /= totFASA
                Ssum /= totFASA
                Lsum /= totFASA
                
                if(Hsum >= 0.4 and Ssum < 0.1):
                    interfaces[int_id]["interface_features"][chain]["secondary_structure_composition"] = "helix"
                elif(Ssum >= 0.4 and Hsum < 0.1):
                    interfaces[int_id]["interface_features"][chain]["secondary_structure_composition"] = "strand"
                elif((Hsum+Ssum) >= 0.4 and Hsum >= 0.1 or Ssum >= 0.1):
                    interfaces[int_id]["interface_features"][chain]["secondary_structure_composition"] = "helix/strand"
                else:
                    interfaces[int_id]["interface_features"][chain]["secondary_structure_composition"] = "irregular"
                interfaces[int_id]["interface_features"][chain]["mean_hydrophobicity_score"] /= interfaces[int_id]["interface_features"][chain]["weight_sum"]
                del interfaces[int_id]["interface_features"][chain]["weight_sum"]
                del interfaces[int_id]["interface_features"][chain]["stack_count"]
                del interfaces[int_id]["interface_features"][chain]["pair_count"]
        
        # Delete non-interacting SSEs from protein data #
        delete = []
        for i in range(len(JSON["protein"]["models"][mi]["secondary_structure_elements"])):
            sse = JSON["protein"]["models"][mi]["secondary_structure_elements"][i]
            interacts = False
            for int_id in interfaces:
                if(sse["id"] in interfaces[int_id]["sse_data"]):
                    interacts = True
                    break
            if(not interacts):
                delete.append(i)
        for i in sorted(delete, reverse=True):
            del JSON["protein"]["models"][mi]["secondary_structure_elements"][i]
        
        # Convert dicts to arrays
        for int_id in interfaces:
            interfaces[int_id]["nucleotide-residue_interactions"] = list(interfaces[int_id]["nucleotide-residue_interactions"].values())
            interfaces[int_id]["nucleotide-sse_interactions"] = list(interfaces[int_id]["nucleotide-sse_interactions"].values())
            interfaces[int_id]["nucleotide_data"] = list(interfaces[int_id]["nucleotide_data"].values())
            interfaces[int_id]["residue_data"] = list(interfaces[int_id]["residue_data"].values())
            interfaces[int_id]["sse_data"] = list(interfaces[int_id]["sse_data"].values())
            interfaces[int_id]["protein_chains"] = list(interfaces[int_id]["protein_chains"])
            interfaces[int_id]["interface_features"] = list(interfaces[int_id]["interface_features"].values())
            for ns in interfaces[int_id]["nucleotide-sse_interactions"]:
                ns["nucleotide_interaction_moieties"] = list(ns["nucleotide_interaction_moieties"])
                ns["residue_interaction_moieties"] = list(ns["residue_interaction_moieties"])
                ns["moiety_interactions"] = list(ns["moiety_interactions"])
        
        # Add list of entity-protein interfaces
        JSON["interfaces"]["models"].append(list(interfaces.values()))
        if(len(interfaces) == 0):
            DELETE_MODELS.append(mi)
    
    ### Delete non-surface residues ###
    for rid in list(JSON["protein"]["residues"].keys()):
        residue = JSON["protein"]["residues"][rid]
        if(not any(residue["surface"])):
            del JSON["protein"]["residues"][rid]
    
    ### Change from dicts to arrays ###
    JSON["protein"]["chains"] = list(JSON["protein"]["chains"].values())
    JSON["protein"]["residues"] = list(JSON["protein"]["residues"].values())
    JSON["dna"]["nucleotides"] = list(JSON["dna"]["nucleotides"].values())
    
    ### Delete any models with no interface ###
    if(len(DELETE_MODELS) == N):
        log("No DNA-protein interactions found.", pdbid)
    mutable_fields = [
        "protein.residues.secondary_structure",
        "protein.residues.sap_score",
        "protein.residues.fasa",
        "protein.chains.secondary_structure",
        "dna.nucleotides.glycosidic_conformation",
        "dna.nucleotides.origin",
        "dna.nucleotides.secondary_structure", 
        "dna.nucleotides.graph_coordinates",
        "protein.models",
        "dna.models",
        "interfaces.models"
    ]
    DELETE_MODELS.sort()
    for mi in DELETE_MODELS:
        mi -= DELETE_MODELS.index(mi)
        deleteModelData(mi, JSON, mutable_fields)
        META_DATA["deleted_models"].append({
            "model": mi,
            "reason": "No DNA-protein interactions were detected."
        })
    
    ### Add meta-data ###
    JSON["meta_data"] = {
        "assembly_chains": META_DATA["assembly_chains"],
        "removed_residues": META_DATA["removed_residues"],
        "added_heavy_atoms": META_DATA["added_heavy_atoms"],
        "deleted_models": META_DATA["deleted_models"],
        "processed_date": str(datetime.date.today()),
        "processing_options": META_DATA["options"],
        "modified_ids": META_DATA["modified_ids"],
        "program_versions": META_DATA["program_versions"],
        "ensemble": META_DATA["ensemble"]
    }
    if(ADD_MMCIF):
        JSON["meta_data"]["citation_data"] = {}
        JSON["meta_data"]["citation_data"]["exp_method"] = mmcif_dict['_exptl.method']
        JSON["meta_data"]["citation_data"]["structure_title"] = mmcif_dict['_struct.title']
        JSON["meta_data"]["citation_data"]["keywords"] = mmcif_dict['_struct_keywords.text']
        if(isinstance(mmcif_dict['_citation.id'], list)):
            index = mmcif_dict['_citation.id'].index('primary')
            JSON["meta_data"]["citation_data"]["citation_title"] = mmcif_dict['_citation.title'][index]
            JSON["meta_data"]["citation_data"]["authors"] = mmcif_dict['_citation_author.name'][index]
            JSON["meta_data"]["citation_data"]["pubmed_id"] = mmcif_dict['_citation.pdbx_database_id_PubMed'][index]
            JSON["meta_data"]["citation_data"]["doi"] = mmcif_dict['_citation.pdbx_database_id_DOI'][index]
            JSON["meta_data"]["citation_data"]["year"] = mmcif_dict['_citation.year'][index]
        else:
            JSON["meta_data"]["citation_data"]["citation_title"] = mmcif_dict['_citation.title']
            JSON["meta_data"]["citation_data"]["authors"] = mmcif_dict['_citation_author.name']
            JSON["meta_data"]["citation_data"]["pubmed_id"] = mmcif_dict['_citation.pdbx_database_id_PubMed']
            JSON["meta_data"]["citation_data"]["doi"] = mmcif_dict['_citation.pdbx_database_id_DOI']
            JSON["meta_data"]["citation_data"]["year"] = mmcif_dict['_citation.year']
        if(isinstance(mmcif_dict['_pdbx_audit_revision_history.revision_date'], list)):
            JSON["meta_data"]["citation_data"]["release_date"] = mmcif_dict['_pdbx_audit_revision_history.revision_date'][0]
        else:
            JSON["meta_data"]["citation_data"]["release_date"] = mmcif_dict['_pdbx_audit_revision_history.revision_date']
    
    ### Write to file ###
    JOUT = open("{}.json".format(pdbid),"w")
    JOUT.write(json.dumps(JSON, indent=None, separators=(',', ':'), sort_keys=True))
    JOUT.close()
    
    ### Add Header info to PDB file ###
    addSSHeader(pdbid, JSON["protein"]["models"][0]["header_info"]["helix_remarks"], JSON["protein"]["models"][0]["header_info"]["sheet_remarks"])
