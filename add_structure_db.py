from pymongo import MongoClient
import json
import os
import sys

def add_structure_db(filepath):
    # MongoDB connection string
    # Adjust the connection string if you have authentication enabled or if you're not running on the default port
    connection_string = "mongodb://localhost:27017/"

    # Path to your JSON file
    # json_file_path = '/home/aricohen/Desktop/dnaprodb/{}.json'.format(filepath)
    json_file_path = filepath

    # Connect to MongoDB
    client = MongoClient(connection_string)

    # Select the database
    db = client['dnaprodb2']

    # Select the collection
    collection = db['dna-protein']

    # Load JSON data from the file
    with open(json_file_path, 'r') as file:
        data = json.load(file)

    # Check if data is a list of documents or a single document
    if isinstance(data, list):
        print("WEIRDWERIDWEIRD\nWEIRDWERIDWEIRD\nWEIRDWERIDWEIRD\nWEIRDWERIDWEIRD\n")
    else:
        # Upsert a single document
        if 'structure_id' not in data:
            print("ERROR. Structure ID not in document for path:", filepath)
            return filepath
        
        structure_id = data["structure_id"]
        if 'external/PDB/pdb_entries/' in structure_id:
            structure_id = structure_id.split('external/PDB/pdb_entries/')[1]
            data['structure_id'] = structure_id

        query = {"structure_id": data["structure_id"]}
        update = {"$set": data}
        result = collection.update_one(query, update, upsert=True)
        if result.matched_count > 0:
            print(f"Updated document with structure_id {data['structure_id']}.")
        elif result.upserted_id is not None:
            print(f"Inserted new document with structure_id {data['structure_id']}.")

    # Close the connection
    client.close()
    return "0"

def get_absolute_filepaths(directory):
    # List to store absolute file paths
    absolute_filepaths = []
    
    # Walk through the directory
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            # Construct absolute path
            absolute_path = os.path.abspath(os.path.join(dirpath, filename))
            absolute_filepaths.append(absolute_path)
    
    return absolute_filepaths

# if __name__ == '__main__':
#     with open('errors.txt', 'w') as my_file:
#         file_paths = get_absolute_filepaths('/home/aricohen/Desktop/dnaprodb/json/')
#         for path in file_paths:
#             output = add_structure_db(path)
#             if not (output == "0"):
#                 my_file.write(output)


if __name__ == '__main__':
    output = add_structure_db(sys.argv[1])