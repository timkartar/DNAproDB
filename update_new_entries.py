#!/usr/bin/env python
import argparse
import os
from os.path import join as ospj
import json
import requests
from string import Template
from pymongo import MongoClient
import gzip
import sys
import shutil

from auto_processStructure import bulkAutoProcessStructures

from dnaprodb_utils import C

def loadConfig(config_file):
    with open(config_file) as FH:
        config = json.load(FH)
    return config

def move_files_to_root(directory, target_directory):
    for root, _, files in os.walk(directory):
        for file in files:
            src_path = os.path.join(root, file)
            dest_path = os.path.join(target_directory, file)
            shutil.move(src_path, dest_path)  # move each file to the target directory

def decompress_gz_files(directory):
    for file_name in os.listdir(directory):
        if file_name.endswith('.gz'):
            file_path = os.path.join(directory, file_name)
            with gzip.open(file_path, 'rb') as f_in:
                with open(file_path[:-3], 'wb') as f_out:  # Remove the last three characters (.gz) for the output file
                    f_out.write(f_in.read())
            os.remove(file_path)  # Optionally remove the .gz file after decompression

# release_date='2023-06-14'
def downloadEntries(release_date):
    CF = loadConfig("annotations_config.json")
    ROOT_DIR = CF["ROOT"]

    ### Query RCSB for new entries
    search_url = CF["RCSB"]["search_url"]
    query = open(ospj(CF["ROOT"], CF["RCSB"]["search_query_template"])).read()
    query = Template(query)
    query = json.loads(query.substitute(release_date='2023-06-14'))

    req = requests.post(search_url, json=query)
    RESULTS = req.json()

    ### Download structures
    root_download_path =ospj(ROOT_DIR, "external/PDB/pdb_entries")
    download_url = CF["RCSB"]["mmcif_data_url"]
    downloaded = open(os.path.join(ROOT_DIR, "external/PDB", "newest_downloaded_releases.txt"), "w")
    for entry in RESULTS["result_set"]:
        entry_id = entry["identifier"].lower()
        path = os.path.join(root_download_path, "{}.cif.gz".format(entry_id))
        # download entry
        div = entry_id[1:3]
        try:
            r = requests.get(download_url.format(div, entry_id))
            with open(path, "wb") as FH:
                FH.write(r.content)
        except Exception as e:
            continue
        downloaded.write(path + '\n')
    downloaded.close()
    decompress_gz_files(root_download_path)
    move_files_to_root(root_download_path, ROOT_DIR)
    return root_download_path

def main(release_date):
    downloadEntries(release_date)
    bulkAutoProcessStructures()

if __name__ == "__main__":
    main('2023-06-14')