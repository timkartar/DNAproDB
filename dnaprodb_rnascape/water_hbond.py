from Bio.PDB import PDBParser, NeighborSearch, Selection, is_aa
import math
import numpy as np

def getWaterHbonds(pdb_filename = 'cifs/2p2r.pdb'):
    # Define criteria for hydrogen bonds
    HBOND_CUTOFF_DISTANCE = 3.5  # Angstroms
    HBOND_ANGLE_CUTOFF = 120     # Degrees
    #COS_ANGLE_CUTOFF = math.cos(math.radians(HBOND_ANGLE_CUTOFF))
    # Load your PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_filename)
    # Separate atoms into protein and water
    protein_atoms = []
    water_atoms = []
    for atom in structure.get_atoms():
        resname = atom.get_parent().get_resname()
        if resname == 'HOH':
            water_atoms.append(atom)
        else:
            protein_atoms.append(atom)
    # Function to calculate the cosine of the angle between two vectors
    def cosine_angle(v1, v2):
        dot_product = sum(a * b for a, b in zip(v1, v2))
        mag_v1 = math.sqrt(sum(a ** 2 for a in v1))
        mag_v2 = math.sqrt(sum(a ** 2 for a in v2))
        return np.arccos(dot_product / (mag_v1 * mag_v2)) * 180/ np.pi
    # Detect water-mediated hydrogen bonds
    hbond_pairs = []
    neighbor_search = NeighborSearch(protein_atoms + water_atoms)
    water_dna = set()
    water_protein = set()
    for water in water_atoms:
        water_coord = water.coord
        nearby_atoms = neighbor_search.search(water_coord, HBOND_CUTOFF_DISTANCE)
        donors = []
        acceptors = []
        for atom in nearby_atoms:
            if atom in protein_atoms:
                vector = atom.coord - water_coord
                if cosine_angle(vector, [0, 0, 1]) >= HBOND_ANGLE_CUTOFF:  # Simplified angle check
                    cand = (water.get_parent().get_id(), (atom.get_parent().get_id(),
                        atom.get_parent().get_parent().id))
                    if is_aa(atom.get_parent()):
                        water_protein.add(cand)
                    else:
                        water_dna.add(cand)
                    '''
                    if atom.element in ['N', 'O']:
                        donors.append(atom)
                    elif atom.element == 'O':
                        acceptors.append(atom)
                    print(donors, acceptors)
                    '''
    water_dna = list(water_dna)
    water_protein  = list(water_protein)
    waters_with_dna = [item[0] for item in water_dna]
    waters_with_protein = [item[0] for item in water_protein]

    waterrr = set(waters_with_dna).intersection(waters_with_protein)

    result = dict()

    for item in waterrr:
        result[item] = [None,None]

    for item in water_dna:
        if item[0] in result.keys():
            result[item[0]][0] = item[1]


    for item in water_protein:
        if item[0] in result.keys():
            result[item[0]][1] = item[1]

    for item in result.values():
        dna = "{}.{}. ".format(item[0][1], item[0][0][1])
        protein = "{}.{}. ".format(item[1][1], item[1][0][1])

        print("{},{}".format(dna, protein))

if __name__ == "__main__":
    getWaterHbonds()
