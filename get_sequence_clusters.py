import os
import json
from pymongo import MongoClient, UpdateOne
from get_current_pdb_ids import get_current_pdb_ids
import requests

resources_link = "."
CLUSTERS = ["30", "40", "50", "70", "90", "95", "100"]
CLUSTER_MAP = {}
PDB_DIR = "/srv/www/dnaprodb.usc.edu/DNAProDB/sequence/clusters/"

def downloadClusters():
    if not os.path.exists(PDB_DIR):
        os.makedirs(PDB_DIR)

    for cluster in CLUSTERS:
        url = f"https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{cluster}.txt"
        local_filename = os.path.join(PDB_DIR, f"clusters-by-entity-{cluster}.txt")

        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(local_filename, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"Downloaded {local_filename}")
        else:
            print(f"Failed to download {url}, status code: {response.status_code}")

def updateExistingDB():
    connection_string = "mongodb://localhost:27017/"
    database_name = "dnaprodb2"
    collection_name = "dna-protein"

    # Connect to MongoDB
    client = MongoClient(connection_string)
    db = client[database_name]
    collection = db[collection_name]

    # Set a batch size for bulk operations
    batch_size = 1000
    operations = []

    for cluster in CLUSTERS:
        path = os.path.join(PDB_DIR, "updateExisting-{}.maps".format(cluster))
        with open(path, 'r') as file:
            for line in file:
                data = json.loads(line)  # Parse the JSON data from each line
                pdbid = data.pop("pdbid")  # Extract and remove the pdbid from the data
                update_query = {"structure_id": pdbid}  # Define the query to find the document
                update_data = {"$set": data}  # Define the update operation

                # Prepare the bulk operation
                operations.append(UpdateOne(update_query, update_data, upsert=True))

                # Execute the batch when it reaches the specified size
                if len(operations) == batch_size:
                    collection.bulk_write(operations)
                    operations = []  # Reset the operations list after executing

    # Perform any remaining operations
    if operations:
        collection.bulk_write(operations)

    client.close()
    print("Updated search metadata for all clusters.")

def updateBcMaps():
    for cluster in CLUSTERS:
        connection_string = "mongodb://localhost:27017/"
        database_name = "dnaprodb2"
        collection_name = "seq{}".format(cluster)  # Replace with your collection name

        # Path to your file
        path = os.path.join(PDB_DIR, "bc-{}.maps".format(cluster))

        # Connect to MongoDB
        client = MongoClient(connection_string)
        db = client[database_name]
        collection = db[collection_name]

        # Ensure indexes are created for quick lookup
        collection.create_index([("pdbid", 1), ("chain_id", 1)])

        operations = []
        # Read the file and insert each line as a document into the collection
        with open(path, "r") as file:
            for line in file:
                try:
                    data = json.loads(line)
                    query = {"pdbid": data["pdbid"], "chain_id": data["chain_id"]}
                    update_data = {"$set": data}
                    operations.append(UpdateOne(query, update_data, upsert=True))
                    # Execute in batches of 1000 (or another suitable number)
                    if len(operations) == 1000:
                        collection.bulk_write(operations)
                        operations = []
                except json.JSONDecodeError as e:
                    print(f"Error decoding JSON: {e}")
                except Exception as e:
                    print(f"An error occurred: {e}")

            try:
                if operations:
                        collection.bulk_write(operations)
            except Exception as e:
                print(f"An error occurred: {e}")
    
        # Close the MongoDB connection
        client.close()
        print("Data insertion complete.")

def createUpdateFile(curUpdates, cluster):
    path = os.path.join(PDB_DIR, "updateExisting-{}.maps".format(cluster))
    with open(path, 'w') as FH:
        for pdbid in curUpdates:
            item = {
                    "pdbid": pdbid.lower(),
                    "search.sequence_clusters.{}".format(cluster): curUpdates[pdbid.lower()],
                }
            FH.write("{}\n".format(json.dumps(item)))

def main():
    downloadClusters()
    for cluster in CLUSTERS:
        with open('{}/clusters-by-entity-{}.txt'.format(PDB_DIR, cluster), 'r') as FH:
            index = 1
            CLUSTER_MAP[cluster] = {}
            for line in FH:
                cid = "{}.{}".format(cluster, index)
                line = line.strip().split()
                for item in line:
                    CLUSTER_MAP[cluster][item.strip()] = cid
                index += 1

    curIds = get_current_pdb_ids()


    # Write cluster mappings to file
    for cluster in CLUSTERS:
        curUpdates = {} # updates to be made if the structure exists in our database
        path = os.path.join(PDB_DIR, "bc-{}.maps".format(cluster))
        with open(path, 'w') as FH:
            for ckey in CLUSTER_MAP[cluster]:
                if (len(ckey) > 10):
                    continue # cover weird edge cases
                pdbid, chain = ckey.split('_')
                cluster_id = CLUSTER_MAP[cluster][ckey]
                item = {
                    "pdbid": pdbid.lower(),
                    "chain_id": chain,
                    "cluster_id": cluster_id,
                    "sequence_identitiy": cluster
                }
                FH.write("{}\n".format(json.dumps(item)))
                if pdbid.lower() in curIds:
                    if pdbid.lower() in curUpdates:
                        curUpdates[pdbid.lower()].append(cluster_id)
                    else:
                        curUpdates[pdbid.lower()] = [cluster_id]
                        
        createUpdateFile(curUpdates, cluster)
        print("file created")

    print("Updating existing DB entries...")
    updateExistingDB()
    print("Updating BC maps into database...")
    updateBcMaps()

if __name__ == '__main__':
    main()
