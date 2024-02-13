from pymongo import MongoClient
import json

# MongoDB connection string
# Adjust the connection string if you have authentication enabled or if you're not running on the default port
connection_string = "mongodb://localhost:27017/"

# Path to your JSON file
json_file_path = '/home/aricohen/Desktop/dnaprodb/2r5z-assembly1.json'

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
    # Upsert multiple documents
    for document in data:
        # Use structure_id as the unique identifier for upsert
        query = {"structure_id": document["structure_id"]}
        update = {"$set": document}
        result = collection.update_one(query, update, upsert=True)
        if result.matched_count > 0:
            print(f"Updated document with structure_id {document['structure_id']}.")
        elif result.upserted_id is not None:
            print(f"Inserted new document with structure_id {document['structure_id']}.")
else:
    # Upsert a single document
    query = {"structure_id": data["structure_id"]}
    update = {"$set": data}
    result = collection.update_one(query, update, upsert=True)
    if result.matched_count > 0:
        print(f"Updated document with structure_id {data['structure_id']}.")
    elif result.upserted_id is not None:
        print(f"Inserted new document with structure_id {data['structure_id']}.")

# Close the connection
client.close()
