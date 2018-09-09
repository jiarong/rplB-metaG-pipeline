# Snakemake workflow: rplB/rpsC-xander-assembly-analyses

[![Snakemake](https://img.shields.io/badge/snakemake-≥4.8.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/jiarong/xander-assembly-pipeline.svg?branch=master)](https://travis-ci.org/jiarong/xander-assembly-pipeline)

This workflow performs a gene targeted ([xander](https://github.com/rdpstaff/Xander_assembler)) assembly of protein coding genes (rplB and rpsC included here) on samples, and generate OTU table and taxonomy table for further microbial diversity analysis.

## Usage

### Step 1: Install workflow (skip if you have `conda` and `snakemake` ready)

This workflow uses conda as package installation tool and snakemake as workflow managment tool. Users just need to install conda and snakemake (via conda), and snakemake will install all dependecies as part of the workflow.

Install conda:

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    conda config --set always_yes yes
    conda update -q conda
    conda info -a
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

Create an conda environment that has snakemake included (in env/xander.yaml):

    conda env create -q --file envs/xander.yaml -n xander
    source activate xander

### Step 2: Configure workflow

    git clone https://github.com/jiarong/xander-assembly-pipeline
    cd xander-assembly-pipelin

Configure the workflow according to your needs via editing the file `config.yaml` and the sheets `metadata.tsv`.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally using `$N` cores:

    snakemake --use-conda --cores $N

Alternatively, it can be run in cluster or cloud environments (see [the docs](http://snakemake.readthedocs.io/en/stable/executable.html) for details).

After successful execution, you will see the OTU table at `PROJECT/output/otu/GENE/otutable.tsv`, and taxonomy table at `PROJECT/output/tax/GENE/taxonomy.tsv` (PRJECT and GENE are defined in `config.yaml`).

