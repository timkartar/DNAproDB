#!/usr/bin/env python
import argparse

from os.path import join as ospj
import json
import requests
from string import Template
from pymongo import MongoClient

from dnaprodb_utils import C

def loadConfig(config_file):
    with open(config_file) as FH:
        config = json.load(FH)
    return config

ROOT_DIR = CF["ROOT"]

### Query RCSB for new entries
search_url = CF["RCSB"]["search_url"]
query = open(ospj(CF["ROOT"], CF["RCSB"]["search_query_template"])).read()
query = Template(query)
query = json.loads(query.substitute(release_date='2019-10-01'))

req = requests.post(search_url, json=query)
RESULTS = req.json()

### Download structures
root_download_path = ospj(ROOT_DIR, "external/PDB/pdb_entries")
download_url = CF["RCSB"]["mmcif_data_url"]
downloaded = open(os.path.join(ROOT_DIR, "external/PDB/pdb_entries", "newest_downloaded_releases.txt"), "w")
for entry in RESULTS["result_set"]:
    entry_id  = entry["identifier"].lower()
    path = os.path.join(root_download_path, entry_id[0], "{}.cif.gz".format(entry_id))
    
    # download entry
    div = entry_id[1:3]
    try:
        r = requests.get(download_url.format(div, entry_id))
        with open(path, "wb") as FH:
            FH.write(r.content)
    except:
        continue
    downloaded.write(path+'\n')
downloaded.close()


# Add data to database
# client = MongoClient()
# db = client[DB_NAME]
# collection = db[COLLECTION_NAME]
