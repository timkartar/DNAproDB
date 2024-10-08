from pymongo import MongoClient
import sys
import json
import requests
from query_jaspar import getJasparLogo
from getUniprot import getUniprot
from get_citation_data import get_citation_data
import time

connection_string = "mongodb://localhost:27017/"

def get_current_pdb_ids():
    client = MongoClient(connection_string)
    db = client['dnaprodb2']
    collection = db['dna-protein']
    
    # Initialize an empty list to hold all the structure_ids
    structure_ids = set()
    
    # Query the collection for all documents and only fetch the structure_id field
    for document in collection.find({}, {"structure_id": 1, "_id": 0}):
        # Check if the structure_id field exists in the document
        if 'structure_id' in document:
            # Add the structure_id to the list
            structure_ids.add(document['structure_id'])
    
    # Close the MongoClient connection
    print(structure_ids)
    return structure_ids

# def delete_all_pdb_ids():
#     # MongoDB connection string
#     connection_string = "mongodb://localhost:27017/"
    
#     # Connect to the MongoDB server
#     client = MongoClient(connection_string)
    
#     # Select the database and collection
#     db = client['dnaprodb2']
#     collection = db['dna-protein']
    
#     # Perform a delete operation where 'pdbid' field exists
#     result = collection.delete_many({"pdbid": {"$exists": True}})
    
#     # Output the result of the deletion
#     print(f"Number of documents deleted: {result.deleted_count}")
    
#     # Close the MongoClient connection
#     client.close()


if __name__ == "__main__":
    get_current_pdb_ids()
