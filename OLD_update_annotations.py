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
from get_citation_data import get_citation_data
import time

connection_string = "mongodb://localhost:27017/"
#NOTE: 1jgg doesn't map correctly

#chid = chain id
def get_entity_id(chid, pdb_id):
    try:
        response = json.loads(requests.get("https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{}/{}".format(pdb_id,
            chid)).text)
        # print(response)
        entity_id = json.dumps(response['rcsb_polymer_entity_instance_container_identifiers']['entity_id'],
                indent=4).replace("\"","")
        return entity_id
    except:
        raise Exception("Could not get entity ID")

def get_uniprot_id(entity_id, pdb_id):
    try:
        response = json.loads(requests.get("https://data.rcsb.org/rest/v1/core/uniprot/{}/{}".format(pdb_id,
            entity_id)).text)
        # print(response)
        uniprot_id = response[0]["rcsb_id"]
        return uniprot_id
    except:
        raise Exception("Could not get Uniprot ID")

def get_all_pdb_ids(client):
    db = client['dnaprodb2']
    collection = db['dna-protein']
    
    # Initialize an empty list to hold all the structure_ids
    structure_ids = []
    
    # Query the collection for all documents and only fetch the structure_id field
    for document in collection.find({}, {"structure_id": 1, "_id": 0}):
        # Check if the structure_id field exists in the document
        if 'structure_id' in document:
            # Add the structure_id to the list
            structure_ids.append(document['structure_id'])
    
    # Close the MongoClient connection
    return structure_ids

# Wrap everything in try to catch errors
def update_annotation(pdb_id, client):
    try:
        # client = MongoClient(connection_string)
        db = client['dnaprodb2']
        collection = db['dna-protein']
        made_change = False

        list_of_ids = []
        id_uniprot_map = {} # maps dnaprodb chain ids to uniprot id

        protein_ids = set()
        uniprot_dict = {}

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
                        protein_ids.add(chain['id'])

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
                entity_id = get_entity_id(id, pdb_id)
                uniprot_id = get_uniprot_id(entity_id, pdb_id)

                # print(id, uniprot_id)
                uniprot_dict[uniprot_id] = True

                # id_uniprot_map[id] = uniprot_id

            except Exception as e:
                continue
        # for search
        protein_names = []
        go_molecular_function_search = []
        organisms = []
        uniprot_ids = []

        for uniprot_id in uniprot_dict:
            # Get Jaspar stuff!
            jaspar_path = ""
            organism = "N/A"
            go_terms_c = []
            go_terms_p = []
            go_terms_f = []
            protein_name = "N/A"
            go_ids=[]
            try:
                jaspar_path = getJasparLogo(uniprot_id)
                if jaspar_path == False:
                    jaspar_path = ""

            except Exception as e:
                print("Error getting jaspar logo!")

            # Get Uniprot data
            try:
                organism, go_terms_c, go_terms_p, go_terms_f, protein_name, go_ids = getUniprot(uniprot_id)
                # print(organism)
                if not go_terms_c:
                    go_terms_c = []
                if not go_terms_f:
                    go_terms_f = []
                if not go_terms_p:
                    go_terms_p = []
            except Exception as e:
                print("Error getting uniprot data")

            uniprot_ids.append({'uniprot_accession': uniprot_id, 'go_ids': list(go_ids)})

            uniprot_object = {}
            uniprot_object['jasparPath'] = jaspar_path
            uniprot_object['organism'] = organism
            uniprot_object['GO_cellular_component'] = [x[0] for x in go_terms_c]
            uniprot_object['GO_biological_process'] =[x[0] for x in go_terms_f]
            uniprot_object['GO_molecular_function'] = [x[0] for x in go_terms_p]
            
            # update search
            for mol_fctn in go_terms_p:
                if mol_fctn not in go_molecular_function_search:
                    go_molecular_function_search.append(mol_fctn)


            uniprot_object['protein_name'] = protein_name
            uniprot_dict[uniprot_id] = uniprot_object
            protein_names.append(protein_name)
            organisms.append(organism)
        search_dict = {}
        # check if uniprot_dict is empty and add dummy stuff
        if not uniprot_dict:
            search_dict['organism'] = '?'
            search_dict['uniprot_names'] = ['?']
            search_dict['GO_molecular_function'] = ['?']
            search_dict['organisms'] = '?'
            search_dict['go_ids'] = ['?']
        else:
            if not organisms:
                organisms = ['?']
            for val in uniprot_dict.values():
                search_dict['organism'] = val['organism']
                break
            if not go_ids:
                search_dict['go_ids'] = ['?']
            search_dict['GO_molecular_function'] = go_molecular_function_search
            search_dict['uniprot_names'] = protein_names
            search_dict['organisms'] = organisms
            search_dict['go_ids'] = list(go_ids)
        search_dict['uniprot_ids'] = uniprot_ids
        document['protein_metadata'] = uniprot_dict
        document['search'] = search_dict

        citation_title, year, authors, doi, pubmed, method, keywords, release_date, title = get_citation_data(pdb_id)
        citation_data = {}
        citation_data['doi'] = doi
        citation_data['structure_title'] = title
        citation_data['release_data'] = release_date
        citation_data['year'] = year
        citation_data['exp_method'] = method
        citation_data['citation_title'] = citation_title
        citation_data['pubmed_id'] = pubmed
        citation_data['authors'] = authors
        citation_data['keywords'] = keywords

        document['meta_data']['citation_data'] = citation_data
        collection.replace_one({"structure_id": pdb_id}, document)
        print(f"Document with structure_id {pdb_id} has been updated.")
    except Exception as e:
        print(f"An error occurred: {e}")

def update_all_annotations():
    try:
        client = MongoClient(connection_string)
        pdb_ids = get_all_pdb_ids(client) # Assume get_all_pdb_ids is modified to accept a client
        for pdb_id in pdb_ids:
            print("Updating", pdb_id)
            update_annotation(str(pdb_id), client)
    except Exception as e:
        print(f"An error occurred during batch update: {e}")
    finally:
        client.close()

def update_one_annotation(pdb_id):
    try:
        client = MongoClient(connection_string)
        print("Updating", pdb_id)
        update_annotation(str(pdb_id), client)
    except Exception as e:
        print(f"An error occurred during update: {e}")
    finally:
        client.close()

if __name__ == "__main__":
    # pdb_id = sys.argv[1]
    # update_annotation(pdb_id)
    update_all_annotations()