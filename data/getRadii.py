import json
import sys
# Elemental vdW radii taken from averages of values in Table 9 of 
# Batsanov S. S., Inorganic Materials, Vol. 37, No. 9, 2001, pp. 871-885.

radii = {
    'element': {
        'LI': (2.2 + 2.63)/2,
        'BE': (1.9 + 2.23)/2,
        'B': (1.8 + 2.05)/2,
        'C': (1.7 + 1.96)/2,
        'N': (1.6 + 1.79)/2,
        'O': (1.55 + 1.71)/2,
        'F': (1.5 + 1.65)/2,
        'NA': (2.4 + 2.77)/2,
        'MG': (2.2 + 2.42)/2, 
        'AL': (2.1 + 2.4)/2,
        'SI': (2.1 + 2.26)/2,
        'P': (1.95 + 2.14)/2,
        'S': (1.8 + 2.06)/2,
        'CL': (1.8 + 2.05)/2,
        'K': (2.8 + 3.02)/2,
        'CA': (2.4 + 2.78)/2,
        'SC': (2.3 + 2.62)/2,
        'TI': (2.15 + 2.44)/2,
        'V': (2.05 + 2.27)/2,
        'CR': (2.05 + 2.23)/2,
        'MN': (2.05 + 2.25)/2,
        'FE': (2.05 + 2.27)/2,
        'CO': (2.0 + 2.25)/2,
        'NI': (2.0 + 2.23)/2,
        'CU': (2.0 + 2.27)/2,
        'ZN': (2.1 + 2.24)/2,
        'GA': (2.1 + 2.41)/2,
        'GE': (2.1 + 2.32)/2,
        'AS': (2.05 + 2.25)/2,
        'SE': (1.9 + 2.18)/2,
        'BR': (1.9 + 2.1)/2,
        'RB': (2.9 + 3.15)/2,
        'SR': (2.55 + 2.94)/2,
        'Y': (2.4 + 2.71)/2,
        'ZR': (2.3 + 2.57)/2,
        'NB': (2.15 + 2.46)/2,
        'MO': (2.1 + 2.39)/2,
        'TC': (2.05 + 2.37)/2,
        'RU': (2.05 + 2.37)/2,
        'RH': (2.0 + 2.32)/2,
        'PD': (2.05 + 2.35)/2,
        'AG': (2.1 + 2.37)/2,
        'CD': (2.2 + 2.37)/2,
        'IN': (2.2 + 2.53)/2,
        'SN': (2.25 + 2.46)/2,
        'SB': (2.2 + 2.41)/2,
        'TE': (2.1 + 2.36)/2,
        'I': (2.1 + 2.22)/2
    }
}

# Additional radii taken from wikipedia https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
radii["element"]["PB"] = 2.02
radii["element"]["HG"] = 1.55

source = sys.argv[1]
count = 0
if(source == 'amber'):
    print("Using AMBER parameters.")
    addToAll = []
    residues = []
    FH = open('amber.siz')
    for line in FH:
        if(count < 9):
            count += 1
            continue
        elif(line[0] == '!' or len(line) < 18):
            count += 1
            continue
        else:
            atm = line[0:5].strip()
            res = line[5:10].strip()
            vdw = float(line[10:18].strip())
            
            if(atm[0].isdigit()):
                continue
            if(atm[0] == 'H'):
                continue
            
            if(len(res) == 0):
                addToAll.append((atm,vdw))
                continue
            
            if(res not in radii):
                residues.append(res)
                radii[res] = {}
            
            radii[res][atm] = vdw
            count += 1
    FH.close()
    for res in residues:
        for add in addToAll:
            radii[res][add[0]] = add[1]
else:
    print("Using NACCESS parameters.")
    FH = open('naccess.config').readlines()
    types = {}
    while count < len(FH):
        line = FH[count]
        if(line.strip() == "types:"):
            # add all defined types enter into sub-loop
            count += 1
            while count < len(FH):
                line = FH[count]
                if(len(line) <= 1):
                    break
                name, radius, _ = line.split()
                types[name] = float(radius)
                count += 1
        else:
            line = line.split('#')[0].split() # strip any comments
            if(len(line) != 3):
                count += 1
                continue
            if(line[0] not in radii):
                radii[line[0]] = {}
            radii[line[0]][line[1]] = types[line[2]]
            count += 1
    # add ANY types to all residues
    ANY = radii.pop("ANY", None)
    for res in radii:
        for atom in ANY:
            radii[res][atom] = ANY[atom]

OUT = open('vdw_radii.json', 'w')
OUT.write(json.dumps(radii, indent=None, separators=(',', ':')))
OUT.close()
