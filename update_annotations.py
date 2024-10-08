from pymongo import MongoClient, UpdateOne
import requests
import json
import logging
from concurrent.futures import ThreadPoolExecutor
import concurrent.futures
from query_jaspar import getJasparLogo
from getUniprot import getUniprot
from get_citation_data import get_citation_data
# NOTE: RUN QUERY_CATH SEPARATELY TO ADD CATH ANNOTATIONS

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def fetch_data(url):
    try:
        response = requests.get(url)
        response.raise_for_status()  # Will raise an exception for HTTP error codes
        return json.loads(response.text)
    except requests.RequestException as e:
        logging.error(f"Failed to fetch data from {url}: {str(e)}")
        return None

def get_entity_id(chid, pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pdb_id}/{chid}"
    response = fetch_data(url)
    if response:
        return response.get('rcsb_polymer_entity_instance_container_identifiers', {}).get('entity_id')
    return None

def get_uniprot_id(entity_id, pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/uniprot/{pdb_id}/{entity_id}"
    response = fetch_data(url)
    if response:
        return response[0].get("rcsb_id")
    return None

def update_annotation(pdb_id, client):
    db = client['dnaprodb2']
    collection = db['dna-protein']
    document = collection.find_one({"structure_id": pdb_id})

    if not document:
        logging.info(f"No entry found for structure_id: {pdb_id}")
        return

    list_of_ids = [chain['id'] for chain in document.get('protein', {}).get('chains', []) if 'id' in chain]
    list_of_ids += [chain['id'] for chain in document.get('dna', {}).get('chains', []) if 'id' in chain]

    uniprot_data = {}
    with ThreadPoolExecutor(max_workers=10) as executor:
        future_to_id = {executor.submit(get_uniprot_id, get_entity_id(id, pdb_id), pdb_id): id for id in list_of_ids}
        for future in concurrent.futures.as_completed(future_to_id):
            id = future_to_id[future]
            uniprot_id = future.result()
            if uniprot_id:
                uniprot_data[uniprot_id] = True  # Placeholder for actual data fetching logic

        protein_names = []
        go_molecular_function_search = []
        organisms = []
        uniprot_ids = []
        genes = []
        try:
            for uniprot_id in uniprot_data:
                # Get Jaspar stuff!
                jaspar_path = ""
                jaspar_id = ["?"]
                organism = "N/A"
                go_terms_c = []
                go_terms_p = []
                go_terms_f = []
                protein_name = "N/A"
                go_ids=[]
                gene_list=[]
                try:
                    jaspar_path, jaspar_id = getJasparLogo(uniprot_id)
                    if jaspar_path == False:
                        jaspar_path = ""
                        jaspar_id = ["?"]
                    else:
                        jaspar_id = [jaspar_id]
                except Exception as e:
                    print("Error getting jaspar logo!")

                # Get Uniprot data
                try:
                    organism, go_terms_c, go_terms_p, go_terms_f, protein_name, go_ids, gene_list = getUniprot(uniprot_id)
                    # print(organism)
                    if not go_terms_c:
                        go_terms_c = []
                    if not go_terms_f:
                        go_terms_f = []
                    if not go_terms_p:
                        go_terms_p = []
                    if not gene_list:
                        gene_list = []
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
                for gene in gene_list:
                    if gene not in genes:
                        genes.append(gene)

                uniprot_object['protein_name'] = protein_name
                uniprot_data[uniprot_id] = uniprot_object
                protein_names.append(protein_name)
                organisms.append(organism)
                
            search_dict = {}
            # check if uniprot_dict is empty and add dummy stuff
            if not uniprot_data:
                search_dict['organism'] = '?'
                search_dict['uniprot_names'] = ['?']
                search_dict['GO_molecular_function'] = ['?']
                search_dict['organisms'] = '?'
                search_dict['go_ids'] = ['?']
                search_dict['genes'] = ['?']
                search_dict['jaspar_id'] = ['?']
            else:
                if not organisms:
                    organisms = ['?']
                for val in uniprot_data.values():
                    search_dict['organism'] = val['organism']
                    break
                if not go_ids:
                    search_dict['go_ids'] = ['?']
                search_dict['GO_molecular_function'] = go_molecular_function_search
                search_dict['uniprot_names'] = protein_names
                search_dict['organisms'] = organisms
                search_dict['go_ids'] = list(go_ids)
                search_dict['genes'] = genes
                search_dict['jaspar_id'] = jaspar_id
            search_dict['uniprot_ids'] = uniprot_ids
            document['protein_metadata'] = uniprot_data
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


    if uniprot_data:
        document['protein_metadata'] = uniprot_data
        collection.replace_one({"structure_id": pdb_id}, document)
        logging.info(f"Document with structure_id {pdb_id} has been updated.")

def update_all_annotations():
    client = MongoClient("mongodb://localhost:27017/")
    db = client['dnaprodb2']
    collection = db['dna-protein']
    structure_ids = collection.distinct("structure_id")

    with ThreadPoolExecutor(max_workers=5) as executor:
        executor.map(lambda pdb_id: update_annotation(pdb_id, client), structure_ids)

    client.close()
    logging.info("Completed updating all annotations.")

def update_single_annotation(pdb_id):
    client = MongoClient("mongodb://localhost:27017/")
    update_annotation(pdb_id, client)

if __name__ == "__main__":
    import sys
    # update_single_annotation(sys.argv[1])
    update_all_annotations()
