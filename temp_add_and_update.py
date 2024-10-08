import sys
import os
from add_structure_db import add_structure_db
from update_annotations import update_single_annotation
from os.path import join as ospj
from pymongo import MongoClient

json_path = "/srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/htdocs/uploads/"

# Connect to MongoDB
client = MongoClient("mongodb://localhost:27017/")
db = client["dnaprodb2"]
collection = db["dna-protein"]

count = 0
for f in os.listdir(json_path):
    pdb = f.split(".")[0]
    # Check if structure_id exists
    if collection.find_one({"structure_id": pdb}):
        print(f"Structure {pdb} already exists in the database. Skipping.")
        continue
    try:
        add_structure_db(ospj(json_path, f))
        update_single_annotation(pdb)
    except Exception as e:
        print("Error occured adding/updating {}\n {}\n".format(pdb, e))
        continue
    count += 1
    print("Addition No.{}".format(count))
    
