from pymongo import MongoClient, UpdateOne
import urllib.request
import os
import gzip
import shutil

save_dir = "./CATH/"
file_name = "cath_domains_list.dat.gz"
# Full path where the file will be saved
save_path = os.path.join(save_dir, file_name)
PDBIDS = {}

def download_mapping():
    # FTP URL, file name to save as, and the directory to save in
    ftp_url = "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/daily-release/newest/cath-b-newest-all.gz"

    # Ensure the save directory exists
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)


    # Download the file from FTP and save it
    try:
        print(f"Downloading {ftp_url} to {save_path}")
        urllib.request.urlretrieve(ftp_url, save_path)
        print("Download complete.")
    except Exception as e:
        print(f"Error downloading file: {e}")
        return

    # Unzip the downloaded file
    try:
        unzipped_file_path = save_path[:-3]  # Remove the '.gz' extension
        with gzip.open(save_path, 'rb') as f_in:
            with open(unzipped_file_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"Unzipped the file to {unzipped_file_path}")
        
        # Optionally, you can remove the .gz file after extraction
        os.remove(save_path)
    except Exception as e:
        print(f"Error unzipping file: {e}")

def process_mapping():
    processed_ids = {}
    print("Reading in CATH data")
    with open('./CATH/cath_domains_list.dat', 'r') as cath_file:
        for line in cath_file:
            id_object = {}
            line = line.split()
            pdbid = line[0][0:4].lower()
            chain = line[0][4]
            cath = line[2].split('.')
            # if(pdbid not in PDBIDS):
            #     continue
            
            Homology = '.'.join(cath[0:4])
            Topology = '.'.join(cath[0:3])
            Architecture = '.'.join(cath[0:2])
            Class = '.'.join(cath[0:1])

            if pdbid in processed_ids:
                if Homology not in processed_ids[pdbid]['Homology']:
                    processed_ids[pdbid]['Homology'].append(Homology)
                if Topology not in processed_ids[pdbid]['Topology']:
                    processed_ids[pdbid]['Topology'].append(Topology)
                if Architecture not in processed_ids[pdbid]['Architecture']:
                    processed_ids[pdbid]['Architecture'].append(Architecture)
                if Class not in processed_ids[pdbid]['Class']:
                    processed_ids[pdbid]['Class'].append(Class)
            else:
                processed_ids[pdbid] = {'Homology': [Homology], 'Topology': [Topology], 'Architecture': [Architecture], 'Class': [Class]}
    return processed_ids
# Call the function
# download_mapping()

def update_database(processed_ids):
    client = MongoClient("mongodb://localhost:27017/")
    db = client['dnaprodb2']
    collection = db['dna-protein']
    
    existing_structure_ids = collection.distinct("structure_id")


    # List to hold all the bulk operations
    bulk_operations = []

    # Loop through each ID in processed_ids
    for pdbid, data in processed_ids.items():
        if pdbid not in existing_structure_ids: # do not insert new structures!
            continue
        # Create an update operation
        update_op = UpdateOne(
            {'structure_id': pdbid},  # Filter to find the document by structure_id
            {
                '$set': {
                    'meta_data.cath': data  # Update or set the 'meta_data.cath' field with the processed data
                }
            },
            upsert=True  # If the document doesn't exist, create it
        )
        # Append the operation to the list
        bulk_operations.append(update_op)

    if bulk_operations:
        # Execute all the operations in bulk
        result = collection.bulk_write(bulk_operations)
        print(f"Bulk update complete. Matched: {result.matched_count}, Modified: {result.modified_count}, Upserted: {result.upserted_count}")
    else:
        print("No operations to perform.")

def main():
    download_mapping()
    processed_ids = process_mapping()
    update_database(processed_ids)
    

if __name__ == "__main__":
    main()
