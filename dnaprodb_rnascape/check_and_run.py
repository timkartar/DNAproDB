import sys
import requests
from xml.etree.ElementTree import fromstring
import json
import wget
import subprocess
import numpy as np

def checkAndRunRNAscape(pdb_id, json_path):
    run = False

    response = requests.get("https://data.rcsb.org/rest/v1/core/entry/{}/".format(pdb_id)).text
    response = json.loads(response)

    assembly_ids = response["rcsb_entry_container_identifiers"]["assembly_ids"]
    entity_ids = response["rcsb_entry_container_identifiers"]["entity_ids"]

    dna_entities = []
    for entity_id in entity_ids:
        response = json.loads(requests.get("https://data.rcsb.org/rest/v1/core/polymer_entity/{}/{}".format(pdb_id,entity_id)).text)
        #print(dict(response).keys())
        try:
            if response["entity_poly"]["rcsb_entity_polymer_type"] == "DNA":
                dna_entities.append(entity_id)
        except:
            continue

    if len(dna_entities) > 2:
        run = True
    else:
        assembly1 = assembly_ids[0]
        response = json.loads(requests.get("https://data.rcsb.org/rest/v1/core/assembly/{}/{}".format(pdb_id,assembly1)).text)
        chains = response["pdbx_struct_assembly_gen"][0]["asym_id_list"]
        ch_count = 0
        for chid in chains:
            try:
                response = json.loads(requests.get("https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{}/{}".format(pdb_id,
                    chid)).text)
                entity_id = json.dumps(response['rcsb_polymer_entity_instance_container_identifiers']['entity_id'],
                       indent=4).replace("\"","")
                if entity_id in dna_entities:
                    ch_count += 1
            except:
                continue
        if ch_count > 2:
            run = True

    if run == True:
        #wget.download("https://files.rcsb.org/download/{}-assembly1.cif".format(pdb_id),
        #        "./cifs/{}.cif".format(pdb_id))
        result = subprocess.run(["python", "run_rnascape.py", "cifs/{}-DNA.pdb".format(pdb_id), pdb_id],
                    capture_output = True,
                    text=True)
        
        l = result.stdout.split("\n")[-2].split(": ")[-1].strip()
        data = dict(np.load(l, allow_pickle=True))
        points = data["points"]
        dssrids = data["dssrids"]

        def convertId(id_string):
            """Converts a DSSR id string to standard DNAproDB format.
            
            Parameters
            ----------
            id_string: string 
                A DSSR nucleotide id string.
            """
            components = id_string.split('.')
            if(components[5] == ''):
                components[5] = ' '

            return '.'.join((components[2], components[4], components[5]))
        
        ids = [convertId(idd) for idd in dssrids]
        j = json.load(open(json_path,"r"))
        
        for i in range(len(j["dna"]["nucleotides"])):
            idd = j["dna"]["nucleotides"][i]["id"]
            idx = ids.index(idd)
            j["dna"]["nucleotides"][i]["graph_coordinates"][0]["circular"]['x'] = points[idx][0]
            j["dna"]["nucleotides"][i]["graph_coordinates"][0]["circular"]['y'] = points[idx][1]

        json.dump(j, open("./updated_{}.json".format(pdb_id),"w"), indent=4)

if __name__ == "__main__":
    pdb_id = sys.argv[1]
    checkAndRunRNAscape(pdb_id, "./4kle.json")
