# DNAproDB Processing Pipeline
The DNAproDB processing pipeline is built as a collection of python scripts that take as input a PDB entry (in mmCIF or PDB file format) and produce a single JSON file as output containing all information about the DNA and protein structural characteristics and the interactions between protein residues and nucleotides in the complex. Additional annotations from several other databases are used as well. A schematic of the pipeline is shown below.

![pipeline](docs/dnaprodb_pipeline.png)

It is preferred to use the mmCIF file format, as the header info is parsed and additonal annotations can be extracted from it, which would be missing for the PDB file format.

## Install Instructions
To install, it is first recommended to create a new python enviroment via conda:

```
conda create -n dnaprodb python=3.8
conda activate dnaprodb
```

Next, clone this repository and install python dependencies via pip

```
git clone https://github.com/jaredsagendorf/dnaprodb-back.git
cd dnaprodb-back
pip install -r python-requirements.txt
```

### Dependencies
In addition to the python packages listed in `python-requirements.txt`, dnaprodb relies on a variety of command line tools that are invoked by the processing scripts. The following must be compiled or made executable and accessible from the system path:

- reduce https://github.com/rlabduke/reduce
- hbplus (source available in `/share/hbplus.tar.gz`)
- Curves 5.3 (source available in `/share/curves5.3_edited.tar.gz`)
- 3DNA (source available in `x3dna-v2.3.tar.gz`)
- pdb2pqr (source available in `pdb2pqr-3.5.2.zip`)
- dssr (available as linux binary in `/share/x3dna-dssr`)
- dssp (available as linux binary in `/share/dssp`)
- SNAP (available as linux binary in `/share/x3dna-snap`)
- msms (available as linux binary in `/share/msms`)

## Run Instructions
To run the pipeline, simply invoke the main script `processStructure.py` from the command line as such:

```
./processStructure.py <FILE_NAME> --clean
```

the output will be a single file containing all the DNAproDB output data in a file called `<FILE_NAME_PREFIX>.json`. 
