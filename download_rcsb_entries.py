#!/usr/bin/env python
import argparse
import os
from os.path import join as ospj
import json
import requests
from string import Template
from pymongo import MongoClient
import datetime
from dnaprodb_utils import C
import auto_processStructure
import query_cath

def loadConfig(config_file):
    with open(config_file) as FH:
        config = json.load(FH)
    return config

CF = loadConfig("annotations_config.json")
ROOT_DIR = CF["ROOT"]

### Query RCSB for new entries
search_url = CF["RCSB"]["search_url"]
query = open(ospj(CF["ROOT"], CF["RCSB"]["search_query_template"])).read()
query = Template(query)
query = json.loads(query.substitute(release_date='1981-03-18'))

req = requests.post(search_url, json=query)
RESULTS = req.json()

### Download structures
root_download_path = ospj(ROOT_DIR, "external/PDB/pdb_entries")
download_url = CF["RCSB"]["mmcif_data_url"]
#downloaded = open(os.path.join(ROOT_DIR, "external/PDB", "newest_downloaded_releases.txt"), "w")

# existing = [f.split(".")[0] for f in os.listdir(ospj(ROOT_DIR, "external/PDB/pdb_entries"))]
existing = set()
with open('ids_considered.txt', 'r') as existing_file:
    for line in existing_file:
        existing.add(line.strip())
entry_id = ""

with open('ids_considered.txt', 'a') as existing_file:
    for entry in RESULTS["result_set"]:
        entry_id  = entry["identifier"].lower()
        print("Considering " + entry_id)
        path = os.path.join('/srv/www/dnaprodb.usc.edu/DNAProDBv3', "{}.cif.gz".format(entry_id))

        # print(path)
        # download entry
        div = entry_id[1:3]
        try:
            if entry_id in existing:
                continue
            print("PROCESSING " + entry_id)
            existing.add(entry_id)
            existing_file.write("\n" + entry_id) 
            auto_processStructure.main(entry_id + ".cif", mPRE_PDB2PQR=True, isUpload=False)
        except Exception as e:
            print(f"Error processing {entry_id}: {str(e)}")
            continue

# add CATH entries
query_cath.main()

today = datetime.date.today()
formatted_date = today.strftime("%B %d, %Y")

with open('last_updated.txt', 'w') as last_updated:
    last_updated.write(formatted_date)
