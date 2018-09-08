import os
import sys
import pandas as pd
from snakemake.utils import validate, min_version

VERSION = '1.0.0'

### set minimum snakemake version ###
min_version("4.8.0")

### load globals

configfile: 'config.yaml'
if config['WORKDIR']:
    workdir: config['WORKDIR']

GENES = config['GENES'].split()
PROJECT = config['PROJECT']
GROUP = config['GROUP']
PATH = config['PATH']
DF = pd.read_table(config['METADATA'])
if GROUP not in DF.columns:
    mes = '*** {} is not a head in {}, please change "GROUP" in config.yaml'
    print(mes.format(GROUP, os.path.basename(config['METADATA'])))
    sys.exit(1)

if PATH not in DF.columns:
    mes = '*** {} is not a head in {}, please change "PATH" in config.yaml'
    print(mes.format(PATH, os.path.basename(config['METADATA'])))
    sys.exit(1)

SAMPLES = DF.loc[:,GROUP].unique().tolist()

config['SAMPLES'] = SAMPLES
config['REF_DIR'] = config['XANDER_DIR']

### check non conda installable dependencies
shell(
     """
     (set +o pipefail; java -jar {config[JAR_DIR]}/hmmgs.jar 2>&1 | grep 'USAGE: HMMgs' > /dev/null || {{ echo "RDPTools are not installed..  Installing.."; wget http://lyorn.idyll.org/~gjr/public2/misc/RDPTools.tgz -O RDPTools.tgz && tar -xzf RDPTools.tgz; }}
        )
     if [[ ! -f Xander_assembler/bin/run_xander_skel.sh ]]; then 
         echo "Xander_assembler not installed.. Installing.."
         (git clone https://github.com/rdpstaff/Xander_assembler.git && cd Xander_assembler/gene_resource && wget https://zenodo.org/record/1410823/files/rpsC.tgz?download=1 -O rpsC.tgz && tar -xzf rpsC.tgz && rm rpsC.tgz)
     fi
     path1=$(echo $PATH | cut -d ':' -f1)
     {config[UCHIME]} 2>&1 | grep "UCHIME 4" > /dev/null || {{ echo "uchime is not installed.. Installing at $path1.."; wget http://drive5.com/uchime/uchime4.2.40_i86linux32 -O $path1/uchime && chmod a+x $path1/uchime; }}
    """
)

# rules
localrules: make_xander_script, make_otutable, make_taxon_table

rule all:
    input:
        #expand('{project}/output/{sample}/k45/{gene}/cluster/{sample}_{gene}_45_final_prot.fasta', project = PROJECT, sample = SAMPLES, gene = GENES),
        #expand('{project}/output/{sample}/k45/{gene}/cluster/{sample}_{gene}_45_coverage.txt', project = PROJECT, sample = SAMPLES, gene = GENES),
        #expand('{project}/output/{sample}/{sample}_xander.sh', project=PROJECT, sample = SAMPLES),
        #expand('{project}/output/{sample}/{sample}_xander_setenv.sh', project = PROJECT, sample = SAMPLES),
        expand('{project}/output/otu/{gene}/otutable.tsv', project = PROJECT, gene = GENES),
        expand('{project}/output/tax/{gene}/taxonomy.tsv', project = PROJECT, gene = GENES),




rule make_xander_script:
    input:
        lambda wildcards: DF.loc[DF.loc[:,GROUP] == wildcards.sample,:].loc[:,PATH],
    output:
        bashfile = expand('{project}/output/{{sample}}/{{sample}}_xander.sh', project = PROJECT)[0],
        envfile = expand('{project}/output/{{sample}}/{{sample}}_xander_setenv.sh', project = PROJECT)[0],
    script:
        "scripts/xander-makebash.py"
        
rule run_xander:
    input:
        bashfile = expand('{project}/output/{{sample}}/{{sample}}_xander.sh', project = PROJECT)[0]
    output:
        alignedprot_list = expand('{project}/output/{{sample}}/k45/{gene}/cluster/{{sample}}_{gene}_45_final_prot_aligned.fasta', project = PROJECT, gene = GENES),
        covfile_list = expand('{project}/output/{{sample}}/k45/{gene}/cluster/{{sample}}_{gene}_45_coverage.txt', project = PROJECT, gene = GENES),
        taxonabunfile_list = expand('{project}/output/{{sample}}/k45/{gene}/cluster/{{sample}}_{gene}_45_taxonabund.txt', project = PROJECT, gene = GENES),
    conda: srcdir('envs/py2.yaml')
    params:
        sourcedir = srcdir('.')
    shell:
        """
        echo "****"
        echo $PATH
        echo "****"
        echo "****"
        echo "****"

        ver=$(python --version 2>&1)
        path1=$(echo $PATH | cut -d ':' -f1)
        [[ $ver == "Python 2"* ]] || echo "*** Xander only work with Python 2; Please add '--use-conda' to snakemake.."
        {config[HMMALIGN]} -h 2>&1 | grep "HMMER 3\." > /dev/null || {{ echo "*** HMMER3 is found for xander.. Please add '--use-conda' to snakemake.."; }}


        for gene in {GENES}; do rm -rf {PROJECT}/output/{wildcards.sample}/k45/${{gene}}; done
        bash {input}
        """
        
rule make_otutable:
    input:
        alignedprot_list = expand('{project}/output/{sample}/k45/{{gene}}/cluster/{sample}_{{gene}}_45_final_prot_aligned.fasta', project = PROJECT, sample = SAMPLES),
        covfile_list = expand('{project}/output/{sample}/k45/{{gene}}/cluster/{sample}_{{gene}}_45_coverage.txt', project = PROJECT, sample = SAMPLES),
    output:
        expand('{project}/output/otu/{{gene}}/otutable.tsv', project = PROJECT)[0]
    script:
        'scripts/make-otutable.py'

rule make_taxon_table:
    input:
        taxonabunfile_list = expand('{project}/output/{sample}/k45/{{gene}}/cluster/{sample}_{{gene}}_45_taxonabund.txt', project = PROJECT, sample = SAMPLES),
    output:
        expand('{project}/output/tax/{{gene}}/taxonomy.tsv', project = PROJECT)[0],
    script:
        'scripts/filter-taxonabund-file.py'



### make percertage identity files
