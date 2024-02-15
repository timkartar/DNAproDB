import requests
import sys
import json

def getUniprot(uniprot):
    response= json.loads(requests.get("https://rest.uniprot.org/uniprotkb/{}?format=json".format(uniprot)).text)
    #print(json.dumps(response, indent=4))

    go_terms_c = [] ## cellular component
    go_terms_f = [] ## biological processes
    go_terms_p = [] ## molecular function
    for item in response["uniProtKBCrossReferences"]:
        database_id = item["database"]
        if database_id == "GO":
            term = item["properties"][0]["value"]
            spl = term.split(":")
            if spl[0] == "C":
                go_terms_c.append(spl[1])
            elif spl[0] == "F":
                go_terms_f.append(spl[1])
            elif spl[0] == "P":
                go_terms_p.append(spl[1])

    organism = response["organism"]["scientificName"]

    return organism, go_terms_c, go_terms_p, go_terms_f

if __name__=="__main__":
    print(getUniprot("Q9RY80"))
