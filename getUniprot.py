import requests
import sys
import json

def getUniprot(uniprot):
    response= json.loads(requests.get("https://rest.uniprot.org/uniprotkb/{}?format=json".format(uniprot)).text)
    #print(json.dumps(response, indent=4))

    go_terms_c = [] ## cellular component
    go_terms_f = [] ## biological processes
    go_terms_p = [] ## molecular function
    go_ids = set()
    for item in response["uniProtKBCrossReferences"]:
        # print(item)
        database_id = item["database"]
        if database_id == "GO":
            go_id = item['id'].split(':')[1]
            term = item["properties"][0]["value"]
            spl = term.split(":")
            if spl[0] == "C":
                go_terms_c.append((spl[1], go_id))
            elif spl[0] == "F":
                go_terms_f.append((spl[1], go_id))
            elif spl[0] == "P":
                go_terms_p.append((spl[1], go_id))
            go_ids.add(go_id)

    organism = response["organism"]["scientificName"]
    try:
        protein_name = response['proteinDescription']['recommendedName']['fullName']['value']
    except:
        print("Could not find protein name!")
        protein_name = "?"

    return organism, go_terms_c, go_terms_p, go_terms_f, protein_name, go_ids

if __name__=="__main__":
    print(getUniprot("P0ACH5"))
