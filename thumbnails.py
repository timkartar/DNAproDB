#!/usr/bin/env python
#updated 6/4/19 12:29 pm
import __main__
__main__.pymol_argv = [ 'pymol' , '-qc'] # Quiet and no GUI
import json
import numpy as np
from Bio.PDB.PDBParser import PDBParser
import sys, time, os, math
import pymol
pymol.finish_launching()
def getCM(model, rids):
    CM = np.array([0,0,0],dtype=float)
    M = 0
    for res_id in rids:
        cid, num, ins = res_id.split('.')
        rid = (' ', int(num), ins)
        if(rid in model[cid]):
            residue = model[cid][rid]
            for atom in residue:
                M += atom.mass
                CM += atom.mass*atom.get_coord()
    CM /= M
    return CM
def getPrincipalAxis(model, rids, CM):
    I = np.zeros((3,3), dtype=float)
    for rid in rids:
        cid, num, ins = rid.split('.')
        rid = (' ', int(num), ins)
        residue = model[cid][rid]
        for atom in residue:
            coord = atom.get_coord() - CM
            mass= atom.mass
            I[0,0] += mass*(coord[1]**2 + coord[2]**2)
            I[0,1] -= mass*coord[0]*coord[1]
            I[0,2] -= mass*coord[0]*coord[2]
            I[1,1] += mass*(coord[0]**2 + coord[2]**2)
            I[1,2] -= mass*coord[1]*coord[2]
            I[2,2] += mass*(coord[0]**2 + coord[1]**2)
    I[1,0] = I[0,1]
    I[2,0] = I[0,2]
    I[2,1] = I[1,2]
    w, v = np.linalg.eig(I)
    return w, v
def fitPlane(model, nids):
    coords = []
    for nid in nids:
        cid, num, ins = nid.split('.')
        nid = (' ', int(num), ins)
        nucleotide = model[cid][nid]
        if('P' in nucleotide):
            coords.append(nucleotide['P'].get_coord())    coords = np.array(coords)
    mean = coords.mean(axis=0)
    coords -= mean
    
    u, s, v = np.linalg.svd(coords)
    n = v[np.argmin(s)]
    
    return n, mean

def getOrthogonal(v1, v2=None):
    if(v2 is None):
        if(v1[2] != 0):
            v = [1, 1, -(v1[0] + v1[1])/v1[2]]
        elif(v1[1] != 0):
            v = [1, -(v1[0] + v1[2])/v1[1], 1]
        else:
            v = [-(v1[1] + v1[2])/v1[0], 1, 1]
        v = np.array(v)/np.linalg.norm(v)
        
        return v
    else:
        v = np.cross(v1, v2)
        v /= np.linalg.norm(v)
        return v
def getAngle(v1, v2):
    return np.arccos(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))*180/np.pi
def chooseAxis(JSON_DATA, model):
    nucleotide_ids = map(lambda x: x["id"], JSON_DATA["dna"]["nucleotides"]) # list of nucleotide ids
    COM = getCM(model, nucleotide_ids) # center of mass of DNA
    num_dna_entities = len(JSON_DATA["dna"]["models"][0]["entities"])
    
    EIG_VALS, EIG_VECS = getPrincipalAxis(model, nucleotide_ids, COM) # principal axis
    ind = np.argsort(EIG_VALS)
    Y = EIG_VECS[ind[1]] # smallest eigenvector
    Z = EIG_VECS[ind[0]]
    X = EIG_VECS[ind[2]]
    
    if(num_dna_entities == 1):
        if(JSON_DATA["dna"]["models"][0]["entities"][0]["type"] in ["perfect_helix", "imperfect_helix", "irregular_helix"]):
            # DNA is helical - use helical axis 
            helix = JSON_DATA["dna"]["models"][0]["entities"][0]["helical_segments"][0]
            if(helix["helical_axis"]["axis_curvature"] == "linear"):
                Y = [
                    helix["helical_axis"]["x_coef"][1],
                    helix["helical_axis"]["y_coef"][1],
                    helix["helical_axis"]["z_coef"][1]
                ]
                Y /= np.linalg.norm(Y)
                X = getOrthogonal(Y)
                Z = getOrthogonal(X, v2=Y)
            else:
                # default to principal axis
                EIG_VALS, EIG_VECS = getPrincipalAxis(model, nucleotide_ids, COM) # principal axis
                ind = np.argsort(EIG_VALS)
                Y = EIG_VECS[ind[1]] # smallest eigenvector
                Z = EIG_VECS[ind[0]]
                X = EIG_VECS[ind[2]]
        else:
            # fit a plane to the DNA phosphates and choose the plane normal
            n, m = fitPlane(model, nucleotide_ids)
    else:
        # For multiple entities, fit a plane to the DNA phosphates and choose the plane normal
        n, m = fitPlane(model, nucleotide_ids)
    
    return X, Y, Z, COM

def rotateVector(vector, axis, angle):
    angle = angle*np.pi/180
    return np.cos(angle)*vector + np.sin(angle)*np.cross(axis, vector) + (1-np.cos(angle))*np.dot(vector, axis)*axis

### Read User Input
STRUCTURE_PATH = os.path.abspath(sys.argv[1])
STRUCTURE_NAME = STRUCTURE_PATH.split('/')[-1].split('.')[0]
DNAPRODB_PATH = os.path.abspath(sys.argv[2])
transparent = "1"
### Load DNAproDB Data File
with open(DNAPRODB_PATH) as FH:
    JSON_DATA = json.load(FH)

### Load Structure in BioPython
bio_parser = PDBParser(QUIET=True)
bio_structure = bio_parser.get_structure("complex", STRUCTURE_PATH)
### Load Structure in Pymol
pymol.cmd.load(STRUCTURE_PATH, STRUCTURE_NAME)
#pymol.cmd.disable("all")
#pymol.cmd.enable(STRUCTURE_NAME)
output_name = sys.argv[1].split('/')[-1].split('.')[0] + "_thumb.png"

#pymol.cmd.set("cartoon_nucleic_acid_color", "orange")
pymol.cmd.show_as("cartoon")
### Determine X-Y-Z axis

# Debug this later...

#X, Y, Z, COM = chooseAxis(JSON_DATA, bio_structure[0])

#### Center the structure around COM
#COM_array = [COM[0], COM[1], COM[2]]
#pymol.cmd.origin(position = COM_array) #COM = [double,double,double]
#pymol.cmd.center("origin")
#### Compute the euler angles between the two coordinate systems
#z = np.array([0.0, 0.0, 1.0])
#y = np.array([0.0, 1.0, 0.0])
#x = np.array([1.0, 0.0, 0.0])#N = np.cross(z, Z)

## First rotation
#c = getAngle(X, N) # gamma (psi)
#pymol.cmd.rotate("z", angle=c, selection='all', camera=0)
#X = rotateVector(X, z, c)
#Y = rotateVector(Y, z, c)
#Z = rotateVector(Z, z, c)

## Second rotation
#b = getAngle(z, Z) # beta (theta)
#pymol.cmd.rotate("x", angle=b, selection='all', camera=0)
#X = rotateVector(X, x, b)
#Y = rotateVector(Y, x, b)
#Z = rotateVector(Z, x, b)

## Third rotation
#N = np.cross(z, Z)
#a = getAngle(x, N) # alpha (phi)
#pymol.cmd.rotate("z", angle=a, selection='all', camera=0)

# Debug this later...
dna_chain_list = []
for chain in JSON_DATA["dna"]["chains"]:
    dna_chain_list.append(chain["id"])
dna_selection = "chain " + "+".join(dna_chain_list)

pymol.cmd.orient(dna_selection)

#change dimensions
pymol.cmd.viewport(240,240)
#turn off reflections and shadows
pymol.cmd.set('ray_shadows','off')
pymol.cmd.set("specular","off")
#pymol.cmd.set("light_count",0)
#turn off depth que
pymol.cmd.set("depth_cue", 0)
###chains must have space. ex: "chain A". ONLY THE PROTEIN chains
pro_chain_list = []
for chain in JSON_DATA["protein"]["chains"]:
    pro_chain_list.append("chain "+chain['id'])

for chain in pro_chain_list:
    pymol.cmd.color("blue", selection = "ss l+'' and " + chain)
pymol.cmd.color("red", selection = "ss H")
pymol.cmd.color("green", selection = "ss S")
pymol.cmd.bg_color("white")
if transparent == "1":
    pymol.cmd.set("ray_opaque_background", 0)
#ensures no clipping... accoring to pymol
pymol.cmd.zoom("all", complete=0)
#pymol.cmd.set("cartoon_nucleic_acid_color", "orange")
#output
pymol.cmd.png(output_name)
pymol.cmd.quit()