# because of temp files, create a new directory and move PDB/CIF file there

# run processStructure.py
# chains > 2, try RNAscape (put this logic in processStructure)
# or, we simply update the coordinates directly in the database!!!!

# if SUCCESS:
# add structure to database
# call update annotations
import processStructure
import sys
import os
from add_structure_db import add_structure_db
from update_annotations import update_single_annotation
import get_sequence_clusters
import json
import shutil
import glob
import traceback

#UPLOAD_PATH = '/home/raktim/dnaprodb.usc.edu/htdocs/uploads'
#UPLOAD_PATH = '/project/rohs_102/raktimmi/dnaprodb.usc.edu/htdocs/uploads'
UPLOAD_PATH = '/srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/htdocs/uploads'

def get_all_file_paths(directory='.'):
    file_paths = []  # List to store file paths
    types = []
    pdbids = []
    for file in os.listdir(directory):
        file_path = os.path.join(directory, file)
        # Check if it's a file and not a directory
        if os.path.isfile(file_path):
            f_split = file.split(".")
            if len(f_split) < 2:
                continue
            type = f_split[1]
            if type in ['cif', 'pdb']:
                file_paths.append(file_path)
                pdbids.append(f_split[0])
                types.append("." + type)
    return pdbids, types, file_paths

def writeFailedStructure(pdbid):
    print("Error found in", pdbid)
    return
    with open("failedStructures.txt", "a") as file:
        file.write(f"{pdbid}\n")
    return "error"

def isFailedStructure(pdbid):
    failedStructures = set()
    with open('./failedStructures.txt', 'r') as file:
        for line in file:
            failedStructures.add(line.strip().lower())
    if pdbid.strip().lower() in failedStructures:
        return True
    return False


def autoProcessStructure(pdbid, type=".pdb", fpath=None, mPRE_PDB2PQR=False, isUpload=False):
    if isFailedStructure(pdbid):
        print("Is failed structure!")
        return "error"
    output = None 
    if fpath is None:
        fpath = "{}{}".format(pdbid, type)
    try:
    	output = processStructure.main(fpath, mPRE_PDB2PQR=mPRE_PDB2PQR)
    except Exception as e:
        print("Exception when processing structure")
        traceback.print_exc()
    json_path = "{}.json".format(pdbid)
    pdb_path = "{}.pdb".format(pdbid)
    if output:
        if os.path.exists(json_path):

            with open(json_path, 'r') as json_file:
                data = json.load(json_file)
                
                if 'error' in data:
                   print("error in data")
                   print(data)
                   writeFailedStructure(pdbid)

                # just move the JSON and PDB files back to uploads, do NOT add to database
                if isUpload:
                    os.chmod(json_path, 0o777)
                    os.chmod(pdb_path, 0o777)    
                    shutil.move(json_path, os.path.join(UPLOAD_PATH, json_path))
                    shutil.move(pdb_path, os.path.join(UPLOAD_PATH, pdb_path))
		    #os.chmod(os.path.join(UPLOAD_PATH, json_path), 777)
                    #os.chmod(os.path.join(UPLOAD_PATH, pdb_path), 777)
                    return
                
            # check if JSON contains error key

            try:
                print("1")
                os.chmod(json_path, 0o777)
                print("2")
                shutil.move(json_path, os.path.join(UPLOAD_PATH, json_path)) ## also keep json
                os.chmod(os.path.join(UPLOAD_PATH, json_path), 0o777)
                print("3 hi Raktim bye")
                if not isUpload:
                    add_structure_db(os.path.join(UPLOAD_PATH, json_path))
                    print("4")
                    update_single_annotation(pdbid)
                    print("5")
            except Exception as e:
                print("7")
                print(e)
                writeFailedStructure(pdbid)
            
        else:
            print("JSON file does not exist for", pdbid)
            writeFailedStructure(pdbid)
    else:
        print(f"Processing failed for {pdbid}")
        writeFailedStructure(pdbid)

def cleanupFiles():
    files_to_delete = ['*.par', '*.dat', 'hstacking.pdb', 'stacking.pdb', 'stacking.json', '*.r3d']

    # Delete specific files
    for pattern in files_to_delete:
        for file in glob.glob(pattern):
            try:
                os.remove(file)
                print(f"Deleted {file}")
            except OSError as e:
                print(f"Error deleting {file}: {e}")

def cleanupAndMove(pdbid, frontendFolder="/srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/htdocs/data", isUpload=True):
    # List of files to delete
    cleanupFiles()

    if not isUpload:
        # Move pdbid.pdb into mvLocation
        mvLocation = os.path.join(frontendFolder, pdbid[-1])
        pdb_file = f"{pdbid}-noH.pdb"

        # Ensure the directory exists
        if not os.path.exists(mvLocation):
            os.makedirs(mvLocation)
            os.system('chmod -R 777 {}'.format(mvLocation))
        # Move the file
        try:
            print(pdb_file)
            shutil.move(pdb_file, os.path.join(mvLocation, pdb_file.replace("-noH","")))
            os.chmod(os.path.join(mvLocation, pdb_file.replace("-noH","")), 0o777)
            print(f"Moved {pdb_file} to {mvLocation}")
        except FileNotFoundError:
            print(f"File {pdb_file} not found. Cannot move.")
        except Exception as e:
            print(f"Error moving {pdb_file}: {e}")
        finally:
            # Delete files starting with pdbid in the current directory
            for file in glob.glob(f"{pdbid}*"):
                try:
                    os.remove(file)
                    print(f"Deleted {file}")
                except OSError as e:
                    print(f"Error deleting {file}: {e}")
    
    # Delete files starting with pdbid in the current directory
    for file in glob.glob(f"{pdbid}*"):
        try:
            os.remove(file)
            print(f"Deleted {file}")
        except OSError as e:
            print(f"Error deleting {file}: {e}")

def bulkAutoProcessStructures():
    cleanupFiles()
    pdbids, types, fpaths = get_all_file_paths()
    for i in range(len(pdbids)):
        print("Running process structure on", pdbids[i], types[1])
        try:
            autoProcessStructure(pdbids[i], types[i])
        except Exception as e:
            print("Error running auto process structure on", pdbids[i], types[1])
        cleanupAndMove(pdbids[i])

def main(filename, mPRE_PDB2PQR, isUpload):
    # NOTE: FIGURE OUT TITLE LATER!!
    pdbid = filename[:-4]

    if not isUpload:
        import subprocess
        subprocess.run(["wget", "https://files.rcsb.org/download/{}".format(filename)])

    file_type = filename[-4:]

    if isUpload:
        autoProcessStructure(pdbid, type=file_type, mPRE_PDB2PQR=mPRE_PDB2PQR, isUpload=True)
        cleanupAndMove(pdbid,
            frontendFolder="/srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/htdocs/uploads/",
            isUpload=True)
    else:
        autoProcessStructure(pdbid, type=file_type, mPRE_PDB2PQR=mPRE_PDB2PQR, isUpload=False)
        cleanupAndMove(pdbid,
            frontendFolder="/srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/htdocs/data/",
            isUpload=False)

if __name__ == '__main__':
    filename = sys.argv[1]
    isUpload = True
    if len(sys.argv) > 3 and sys.argv[3].strip() == 'update':
        isUpload = False
    mPRE_PDB2PQR = bool(int(sys.argv[2]))
    main(filename, mPRE_PDB2PQR, isUpload)
