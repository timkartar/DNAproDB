from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBParser, MMCIFParser
import sys

parser=MMCIFParser()
f = sys.argv[1]
model = parser.get_structure("2r5z",f)[0]
dssp = DSSP(model, f)
print(dict(dssp))
