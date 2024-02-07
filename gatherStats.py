#!/usr/bin/env python

import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.4f')
import sys
import os
import glob
import copy
import math
import numpy as np
from dnaprodb_utils import C
#import matplotlib.pyplot as plt
DATA_PATH = C["DATA_PATH"]

# Load PDB Components dictionary
with open(os.path.join(DATA_PATH,'components.json')) as FILE:
    COMPONENTS = json.load(FILE)

def sumDicts(acc, val, scale):
    for key in acc:
        if(key not in val):
            continue
        if(isinstance(acc[key], dict)):
            sumDicts(acc[key], val[key], scale)
        else:
            acc[key] += val[key]/scale

def makeStats(data, nonNumericKeys, nonNumericValues, key_string=[]):
    for key in data:
        if(isinstance(data[key], dict)):
            makeStats(data[key], nonNumericKeys, nonNumericValues, key_string+[key])
        elif(isinstance(data[key], list)):
            if(key in nonNumericKeys):
                # get counts of array of values
                values = nonNumericValues[nonNumericKeys.index(key)]
                counts = {}
                for value in values:
                    if(len(data[key]) > 0):
                        counts[value] = float(data[key].count(value))/len(data[key])
                    else:
                        counts[value] = 0.0
                data[key] = counts
            else:
                # get stats of array of numbers
                m = np.array(data[key], dtype=np.float32)
                cutoff_lower = 0.5
                cutoff_upper = 999.9
                ks = '.'.join(key_string+[key])
                
                #bw = 0.3
                #sb = cutoff_lower
                #a = a[np.logical_or(a >= cutoff_lower, a < cutoff_upper)]
                #nbins = 0
                
                #if( len(a) > 0 ):
                    ## Determine a proper cut-off for masking: Want to remove upper and lower outlier bins
                    #nbins = int(math.ceil((a.max()-a.min())/bw))
                    #if(nbins > 9):
                        #h, e = np.histogram(a, nbins)
                        #h = h[h > 0]
                        #mh = np.median(h)
                        #ma = np.median(a)
                        #MAD = np.median(np.abs(h-mh))
                        #if(MAD > 0):
                            #Z = np.abs(0.6745*(h-mh)/MAD)
                            ## Set bottom tail cut-off
                            #iO = -1 # index of previously seen outlier
                            #for i in range(len(Z)):
                                #if(Z[i] > 3.5):
                                    #if(i-iO > 1):
                                        #break
                                    #cutoff_lower = (i+1)*bw + sb
                            
                            ## Set upper tail cut-off
                            #iO = len(Z)+1
                            #for i in reversed(range(len(Z))):
                                #if(Z[i] > 3.5):
                                    #if(iO- i > 1):
                                        #break
                                    #cutoff_upper = (i-1)*bw + sb
                
                ## Mask the values removing outliers
                #print(nbins, len(a), cutoff_lower, cutoff_upper, ks)
                #if(cutoff_lower == cutoff_upper):
                    #print(h)
                    #print(Z)
                    #print(mh, ma, MAD)
                    #plt.hist(a, nbins, normed=1, facecolor='g', alpha=0.75)
                    #plt.title(ks)
                    #plt.show()
                    #exit(0)
                
                m = m[np.logical_and(m >= cutoff_lower, m < cutoff_upper)]
                if(m.size > 0):
                    stats = {
                        "mean": m.mean(),
                        "std": m.std(),
                        "median": np.median(m),
                        "max": m.max(),
                        "min": m.min(),
                        "p10": np.percentile(m, 10),
                        "p20": np.percentile(m, 20),
                        "p50": np.percentile(m, 50),
                        "p80": np.percentile(m, 80),
                        "p90": np.percentile(m, 90)
                    }
                    #if(ks in ("DGARG.basa.sr.sc", "DALYS.basa.pp.sc", "DCHIS.vdw_sum.pp.sc", "DTTYR.hbond_sum.sr.sc")):
                        #print(m.min())
                        #print(m.max())
                        #bw = 0.5
                        #nbins = int(math.ceil((m.max()-m.min())/bw))
                        #if(nbins > 1):
                            #plt.hist(m, nbins, normed=1, facecolor='g', alpha=0.75)
                            #plt.title(ks)
                            #plt.show()
                else:
                    stats = {
                        "mean": 0,
                        "std": 0,
                        "median": 0,
                        "max": 0,
                        "min": 0,
                        "p10": 0,
                        "p20": 0,
                        "p50": 0,
                        "p80": 0,
                        "p90": 0
                    }
                for s in stats:
                    stats[s] = float(stats[s])
                data[key] = stats
        else:
            continue

nucMtyLabel = ['wg', 'sg', 'bs', 'sr', 'pp']
resMtyLabel = ['mc', 'sc']

# Interactions
int_mty_fields = ["basa", "hbond_sum", "vdw_sum"]
int_template = {
    "geometry": [],
    "mean_nn_distance": [],
    "min_distance": [],
    "cm_distance": [],
    "count": 0
}
for f in int_mty_fields:
    int_template[f] = {}
    for nmty in nucMtyLabel:
        int_template[f][nmty] = {}
        for rmty in resMtyLabel:
            int_template[f][nmty][rmty] = []

# Nucleotides
nuc_mty_fields = ["hbond_sum", "vdw_interaction_sum"]
nuc_template = {
    "basa_sum": {},
    "residue_interaction_count": [],
    "count": 0
}
for nmty in nucMtyLabel:
    nuc_template["basa_sum"][nmty] = []
for f in nuc_mty_fields:
    nuc_template[f] = {}
    for nmty in nucMtyLabel:
        nuc_template[f][nmty] = {}
        for rmty in resMtyLabel:
            nuc_template[f][nmty][rmty] = []

# Residues
res_mty_fields = ["hbond_sum", "vdw_interaction_sum"]
res_template = {
    "basa_sum": {},
    "nucleotide_interaction_count": [],
    "count": 0
}
for rmty in resMtyLabel:
    res_template["basa_sum"][rmty] = []
for f in res_mty_fields:
    res_template[f] = {}
    for nmty in nucMtyLabel:
        res_template[f][nmty] = {}
        for rmty in resMtyLabel:
            res_template[f][nmty][rmty] = []

# Data Dict
nucleotides = ['DA','DC','DG','DT','DI','DU']
residues = [
    'ALA','ARG','ASN','ASP','CYS',
    'GLU','GLN','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL'
]
DATA = {} # keys: NUC-RES, RES, NUC
for nuc in nucleotides:
    DATA[nuc] = copy.deepcopy(nuc_template)
    for res in residues:
        DATA[nuc+res] = copy.deepcopy(int_template)
for res in residues:
    DATA[res] = copy.deepcopy(res_template)

# loop over directories
dirlist = sys.argv[1]
DIRS = open(dirlist)
for d in DIRS:
    d = d.strip()
    for listFile in glob.glob("./{}/list*.txt".format(d)):
        LF = open(listFile)
        for pdbid in LF:
            pdbid = pdbid.strip()
            if(os.access("{}/{}.json".format(d,pdbid), os.R_OK)):
                with open("{}/{}.json".format(d,pdbid)) as D:
                    data = json.load(D)
                if("error" in data):
                    continue
                
                # Make nuc and res lookups
                nucDict = {}
                resDict = {}
                try:
                    for n in data["dna"]["nucleotides"]:
                        nucDict[n["id"]] = n
                    for r in data["protein"]["residues"]:
                        resDict[r["id"]] = r
                except:
                    print(("{}: json data has unexpected format.".format(pdbid)))
                    continue
                # Loop over model data
                for model in data["interfaces"]["models"]:
                    for interface in model:
                        for nr in interface["nucleotide-residue_interactions"]:
                            res = nr["res_name"]
                            if(res not in residues):
                                res = COMPONENTS[res]['_chem_comp.mon_nstd_parent_comp_id']
                            nuc = nr["nuc_name"]
                            if(nuc not in nucleotides):
                                nuc = COMPONENTS[nuc]['_chem_comp.mon_nstd_parent_comp_id']
                            
                            # Add moiety fields
                            for f in int_mty_fields:
                                for nmty in nucMtyLabel:
                                    if(not (nmty in nr[f])):
                                        continue
                                    for rmty in resMtyLabel:
                                        DATA[nuc+res][f][nmty][rmty].append(nr[f][nmty][rmty])
                            
                            # Add geometry
                            DATA[nuc+res]["geometry"].append(nr["geometry"])
                            
                            # Add min_dist
                            DATA[nuc+res]["min_distance"].append(nr["min_distance"])
                            
                            # Add mean_nn_dist
                            DATA[nuc+res]["mean_nn_distance"].append(nr["mean_nn_distance"])
                            
                            # Add cm distance
                            DATA[nuc+res]["cm_distance"].append(nr["cm_distance"])
                            
                            # Add count
                            DATA[nuc+res]["count"] += 1
                        
                        for n in interface["nucleotide_data"]:
                            nuc = nucDict[n["nuc_id"]]["name"]
                            if(nuc not in nucleotides):
                                nuc = COMPONENTS[nuc]['_chem_comp.mon_nstd_parent_comp_id']
                            
                            # Add moiety fields
                            for f in nuc_mty_fields:
                                for nmty in nucMtyLabel:
                                    if(not (nmty in n[f])):
                                        continue
                                    for rmty in resMtyLabel:
                                        DATA[nuc][f][nmty][rmty].append(n[f][nmty][rmty])
                            
                            # Add BASA
                            for nmty in nucMtyLabel:
                                if(not (nmty in n["basa_sum"])):
                                    continue
                                DATA[nuc]["basa_sum"][nmty].append(n["basa_sum"][nmty])
                            
                            # Add interaction counts
                            DATA[nuc]["residue_interaction_count"].append(n["residue_interaction_count"])
                            DATA[nuc]["count"] += 1
                        
                        for r in interface["residue_data"]:
                            res = resDict[r["res_id"]]["name"]
                            if(res not in residues):
                                res = COMPONENTS[res]['_chem_comp.mon_nstd_parent_comp_id']
                            
                            # Add moiety fields
                            for f in res_mty_fields:
                                for nmty in nucMtyLabel:
                                    if(not (nmty in r[f])):
                                        continue
                                    for rmty in resMtyLabel:
                                        DATA[res][f][nmty][rmty].append(r[f][nmty][rmty])
                            
                            # Add BASA
                            for rmty in resMtyLabel:
                                if(not (rmty in r["basa_sum"])):
                                    continue
                                DATA[res]["basa_sum"][rmty].append(r["basa_sum"][rmty])
                            
                            # Add interaction counts
                            DATA[res]["nucleotide_interaction_count"].append(r["nucleotide_interaction_count"])
                            DATA[res]["count"] += 1
        LF.close()
DIRS.close()

# Compile statistics
nonNumericKeys = ["geometry"]
nonNumericValues = [["psuedo_pair", "psuedo_stack", "none"]]
makeStats(DATA, nonNumericKeys, nonNumericValues)

# Add fallback stats
DATA["fallback"] = {}
for f in int_mty_fields:
    DATA["fallback"][f] = {}
    for n in nucMtyLabel:
        DATA["fallback"][f][n] = {}
        for r in resMtyLabel:
            DATA["fallback"][f][n][r] = {
                "mean": 0,
                "std": 0,
                "median": 0,
                "max": 0,
                "min": 0,
                "p10": 0,
                "p20": 0,
                "p50": 0,
                "p80": 0,
                "p90": 0
            }

# Add stat averages
for n in nucleotides[:-2]:
    for r in residues:
        sumDicts(DATA["fallback"], DATA[n+r], 80.0)

OUT = open(os.path.join(DATA_PATH, "interaction_stats.json"), "w")
OUT.write(json.dumps(DATA, indent=None,separators=(',', ':'), sort_keys=True))
OUT.close()
