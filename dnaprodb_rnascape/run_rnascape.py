# Quick Example: 'python run_rnascape.py 1ivs.cif 1ivs'
# This script runs RNAscape on a cif or pdb file using your specified settings
# If downloading from the PDB, use a biological assembly file

# BEGIN USER SETTINGS ----------------------------------------------
# NOTE: file must be cif or pdb
import sys
import os

cif = sys.argv[1].strip()
prefix = sys.argv[2].strip()
X3DNA_PATH = 'x3dna-dssr' # path to DSSR binary (e.g., /bin/x3dna-dssr')
json_path = os.path.dirname(os.path.abspath(__file__)) + '/json/' # by default, DSSR output will be saved to the json folder
MEDIA_PATH = os.path.dirname(os.path.abspath(__file__))+'/output/' # where all output of RNAscape will be stored
FIG_PATH = '/processed_images/' # figures will be in MEDIA_PATH/FIG_PATH

# default arguments for RNAscape
cond_bulging = False # this bulges out all loops to create more space on the graph. Change to True to only conditionally bulge loops (condense them)

bp_type = "dssrLw" #BP annotation types: dssrLW, dssr, saenger, none
# dssrLw refers to Leontis-Westhof base pair annotation style

rotation = 0 # float between 0 and 360
extra={'arrowsize':1, 'circlesize':1,
            'circle_labelsize':1, 'cols':['#FF9896', '#DBDB8D', '#90CC84', '#AEC7E8', '#FFFFFF'], #A,C,G,U,X (non-standard)
            'showNumberLabels': False, 'numberSeparation': 1, 'numberSize': 1, 'counter': 0, 'markerSize': 1,
            'numberColor': 'saddlebrown', 
} # other small visual settings
# END USER SETTINGS -------------------------------------------------

from rnascape import rnascape
import os, sys, subprocess
from plot import Plot
import time
import random
import numpy as np

# run DSSR and save output to json folder
subprocess.run([X3DNA_PATH,"-i={}".format(cif),"-o={}/{}-dssr.json".format(json_path, prefix),"-idstr=long","--json","--prefix={}".format(prefix)], cwd=os.getcwd())
subprocess.run([X3DNA_PATH,"-i={}".format(cif),"-o={}/{}-dssr.json".format(json_path, prefix),"-idstr=long","--json","--prefix={}".format(prefix),
    "--cleanup"], cwd=os.getcwd())
points, markers, ids, chids, dssrids, dssrout, prefix = rnascape(prefix, cif, '{}/{}-dssr.json'.format(json_path, prefix),
        cond_bulging=cond_bulging, mFIG_PATH=FIG_PATH, mDSSR_PATH = X3DNA_PATH)

# to keep track of each file
time_string = str(int(time.time())) + str(random.randint(0,100))

figpath, pngpath, log = Plot(points, markers, ids, chids, dssrids, dssrout, prefix, 
                             bp_type=bp_type, extra=extra, time_string=time_string, 
                             rotation=rotation, mFIG_PATH=FIG_PATH, mMEDIA_PATH=MEDIA_PATH)

# Uncomment below lines to save geometrically mapped points and other output to an NPZ file
# npz_filepath = "{}/saved_output/{}.npz".format(MEDIA_PATH,time_string)
# np.savez(npz_filepath, points=points, markers=markers, ids=ids, dssrids=dssrids,
#         chids=chids, prefix=prefix, bp_type=bp_type, extra=extra,
#         time_string=time_string, log=log)

with open("{}/saved_output/{}.log".format(MEDIA_PATH,time_string), 'w') as log_file:
    log_file.write(log)

# print non-WC nucleotides if applicable
print(log)
print('Output:', MEDIA_PATH + figpath)
print('Output:', MEDIA_PATH + pngpath)

# Save output of rnascape function to enable regeneration of labels!
npz_filepath = "{}/saved_output/{}.npz".format(MEDIA_PATH,time_string)

np.savez(npz_filepath, points=points, markers=markers, ids=ids, dssrids=dssrids,
         chids=chids, prefix=prefix, bp_type=bp_type, extra=extra,
         time_string=time_string, log=log)

print("saved raw outputs at: " + npz_filepath)

