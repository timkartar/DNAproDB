import sys
import os
from add_structure_db import add_structure_db
from update_annotations import update_single_annotation
from os.path import join as ospj

json_path = "/srv/www/dnaprodb.usc.edu/data2024_2/uploads/"

old_json_path = "/srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/htdocs/uploads/"
count = 0
for f in os.listdir(json_path):
    if not os.path.exists(ospj(old_json_path,f)):
        pdb = f.split(".")[0]
        try:
            add_structure_db(ospj(json_path, f))
            update_single_annotation(pdb)
        except Exception as e:
            print("Error occured adding/updating {}\n {}\n".format(pdb, e))
            continue
        count += 1
        print("Addition No.{}".format(count))
    
