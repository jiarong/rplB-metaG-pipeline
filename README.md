# Snakemake workflow: rplB/rpsC-xander-assembly-analyses

[![Snakemake](https://img.shields.io/badge/snakemake-≥4.8.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/rna-seq-star-deseq2.svg?branch=master)](https://travis-ci.org/snakemake-workflows/rna-seq-star-deseq2)

This workflow performs a gene targeted (xander) assembly of rplB on samples and generate OTU table and taxonomy table.

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/snakemake-workflows/rna-seq-spew/releases).
If you intend to modify and further develop this workflow, fork this reposity. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, once available, its DOI.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml` and the sheets `metadata.tsv`.

### Step 3: Execute workflow

All you need to execute this workflow is to install Snakemake via the [Conda package manager](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda). Software needed by this workflow is automatically deployed into isolated environments by Snakemake.

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores. Alternatively, it can be run in cluster or cloud environments (see [the docs](http://snakemake.readthedocs.io/en/stable/executable.html) for details).

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.

### Step 4: Investigate results

After successful execution, you can create a self-contained report with all results via:

    snakemake --report report.html
