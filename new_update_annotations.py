#FLOW:
# Given PDB ID:
# Open Database and pull structure
# Get Chain IDs
    # Get Entity ID
        # Get Uniprot ID
            #From here call a lot of code!
            # Jaspar
            # etc.
# Get citation data
# Add date last modified


from pymongo import MongoClient
import sys
import json
import requests
from query_jaspar import getJasparLogo
from getUniprot import getUniprot

connection_string = "mongodb://localhost:27017/"
pdb_id = '3ldy'
#NOTE: 1jgg doesn't map correctly

#chid = chain id
def get_entity_id(chid):
    try:
        response = json.loads(requests.get("https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{}/{}".format(pdb_id,
            chid)).text)
        # print(response)
        entity_id = json.dumps(response['rcsb_polymer_entity_instance_container_identifiers']['entity_id'],
                indent=4).replace("\"","")
        return entity_id
    except:
        raise Exception("Could not get entity ID")

def get_uniprot_id(entity_id):
    try:
        response = json.loads(requests.get("https://data.rcsb.org/rest/v1/core/uniprot/{}/{}".format(pdb_id,
            entity_id)).text)
        # print(response)
        uniprot_id = response[0]["rcsb_id"]
        return uniprot_id
    except:
        raise Exception("Could not get Uniprot ID")

try:
    client = MongoClient(connection_string)
    db = client['dnaprodb2']
    collection = db['dna-protein']

    list_of_ids = []
    id_uniprot_map = {} # maps dnaprodb chain ids to uniprot id

    # Replace 'your_structure_id' with the actual structure_id you're looking for
    document = collection.find_one({"structure_id": pdb_id})

    if document is None:
        print(f"No entry found for structure_id: {pdb_id}")
        sys.exit(0)  # Exit gracefully if the document is not found

    if 'protein' in document:
        protein = document['protein']
        if 'chains' in protein:
            chains = protein['chains']
            for chain in chains:
                if 'id' in chain:
                    list_of_ids.append(chain['id'])
        else:
            print("Error: 'chains' list not found in the 'protein' object.")
    else:
        print("Error: 'protein' object not found in the document.")

    if 'dna' in document:
        dna = document['dna']
        if 'chains' in dna:
            chains = dna['chains']
            for chain in chains:
                if 'id' in chain:
                    list_of_ids.append(chain['id'])
        else:
            print("Error: 'chains' list not found in the 'dna' object.")
    else:
        print("Error: 'dna' object not found in the document.")

    for id in list_of_ids:
        try:
            entity_id = get_entity_id(id)
            uniprot_id = get_uniprot_id(entity_id)
            print(id, uniprot_id)
            id_uniprot_map[id] = uniprot_id
        except:
            continue
    
    for id, uniprot_id in id_uniprot_map.items():
        
        # Get Jaspar stuff!
        jaspar_path = None
        organism = None
        go_terms_c = None
        go_terms_p = None
        go_terms_f = None
        try:
            jaspar_path = getJasparLogo(uniprot_id)
            if jaspar_path == False:
                jaspar_path = None
        except:
            print("Error getting jaspar logo!")

        # Get Uniprot data
        try:
            organism, go_terms_c, go_terms_p, go_terms_f = getUniprot(uniprot_id)
            print(organism)
            if not go_terms_c:
                go_terms_c = None
            if not go_terms_f:
                go_terms_f = None
            if not go_terms_p:
                go_terms_p = None
        except:
            print("Error getting uniprot data")



except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit(1)  # Exit gracefully on any other error
