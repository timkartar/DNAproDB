def getJasparLogo(uniprot):
    lines = [l.strip() for l in open("./JASPAR/uniprot_to_jaspar_map.tsv","r").readlines()]
    mapd = {}
    for item in lines:
        spl = item.split("\t")
        mapd[spl[0]] = spl[1]
    if uniprot in mapd.keys():
        return "/JASPAR/logos/{}.svg".format(mapd[uniprot])
    else:
        return False
if __name__=="__main__":
    print(getJasparLogo('P17676'))
