#!/usr/bin/bash
export PATH="/srv/www/dnaprodb.usc.edu/DeepPBS/dependencies/bin:$PATH"
export X3DNA='/srv/www/dnaprodb.usc.edu/DeepPBS/x3dna-v2.3-linux-64bit/x3dna-v2.3'

export PATH="/srv/www/dnaprodb.usc.edu/DNAProDB/share:$PATH"
export PATH="/srv/www/dnaprodb.usc.edu/DNAProDB/share/hbplus:$PATH"

cd /srv/www/dnaprodb.usc.edu/DNAProDB
nohup /srv/www/dnaprodb.usc.edu/miniconda3/envs/dnaprodb/bin/python download_rcsb_entries.py > update_log.txt
mongodump --db dnaprodb2 --out /srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/htdocs/data
tar -czvf /srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/htdocs/data/dnaprodb2.tar.gz /srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/htdocs/data/dnaprodb2

# Directory to search for files
directory="/srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/htdocs/uploads"

# Find and remove files with names 10 characters or longer (excluding extension)
find "$directory" -type f | awk -F/ '{name=$(NF); split(name, arr, /\./); if (length(arr[1]) >= 10) print}' | xargs rm -f

touch last_run_date.txt
