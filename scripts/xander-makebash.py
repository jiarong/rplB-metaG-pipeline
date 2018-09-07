import sys
import os
import glob
import errno

GENES = snakemake.config['GENES']
#
# Parameters to change in xander_setenv.sh
#
REF_DIR = os.path.abspath(os.path.expanduser(snakemake.config['REF_DIR']))
JAR_DIR = os.path.abspath(os.path.expanduser(snakemake.config['JAR_DIR']))
UCHIME = snakemake.config['UCHIME']
HMMALIGN = snakemake.config['HMMALIGN']

## THIS SECTION MUST BE MODIFIED BASED ON THE INPUT DATASETS
## De Bruijn Graph Build Parameters
K_SIZE = snakemake.config['K_SIZE']  # kmer size, should be multiple of 3 = 4 GB, increase FILTER_SIZE if the bloom filter predicted false positive rate i s greater than 1%
MIN_COUNT = snakemake.config['MIN_COUNT']  # minimum kmer abundance in SEQFILE to be included in the final de Bruijn graph structure
THREADS = snakemake.config['THREADS']

FILTER_SIZE = snakemake.config['FILTER_SIZE']  # memory = 2**(FILTER_SIZE-2), 38 = 64 GB, 37 = 32 GB, 36 = 16 GB
MAX_JVM_HEAP = snakemake.config['MAX_JVM_HEAP']  # memory for java program, must be larger than the corresponding memory of the FILTER_SIZE

def main():
    f_list = [os.path.abspath(f) for f in snakemake.input]
    print(f_list)
    group = snakemake.wildcards.sample
    envfile = snakemake.output['envfile']
    qsub = snakemake.output['bashfile']

    wkdir = os.path.abspath(os.path.dirname(envfile))
    #try:
    #    os.mkdir(wkdir)
    #except OSError as e:
    #    if e.errno != errno.EEXIST:
    #        raise
    template = os.path.join(REF_DIR, 'bin', 'xander_setenv.sh')
    with open (envfile, 'w') as fw_envfile:
        for line in open(template):
            if line.startswith('SEQFILE'):
                line = 'SEQFILE="{}"\n'.format(' '.join(f_list))
            if line.startswith('WORKDIR'):
                line = 'WORKDIR="{}"\n'.format(wkdir)
            if line.startswith('REF_DIR'):
                line = 'REF_DIR="{}"\n'.format(REF_DIR)
            if line.startswith('JAR_DIR'):
                line = 'JAR_DIR="{}"\n'.format(JAR_DIR)
            if line.startswith('UCHIME'):
                line = 'UCHIME="{}"\n'.format(UCHIME)
            if line.startswith('HMMALIGN'):
                line = 'HMMALIGN="{}"\n'.format(HMMALIGN)
            if line.startswith('SAMPLE_SHORTNAME'):
                line = 'SAMPLE_SHORTNAME="{}"\n'.format(group)
            if line.startswith('K_SIZE'):
                line = 'K_SIZE="{}"\n'.format(K_SIZE)
            if line.startswith('FILTER_SIZE'):
                line = 'FILTER_SIZE="{}"\n'.format(FILTER_SIZE)
            if line.startswith('MAX_JVM_HEAP'):
                line = 'MAX_JVM_HEAP="{}"\n'.format(MAX_JVM_HEAP)
            if line.startswith('MIN_COUNT'):
                line = 'MIN_COUNT="{}"\n'.format(MIN_COUNT)

            fw_envfile.write(line)


    with open(qsub, 'w') as fw_qsub:
        print('Xanderdir={}'.format(REF_DIR), file=fw_qsub)
        print('Setenvfile={}'.format(envfile), file=fw_qsub)
        print('Genes="{}"'.format(GENES), file=fw_qsub)
        
        print('set -e', file=fw_qsub)
        
        print('echo *** "build starting.."', file=fw_qsub)
        print('bash $Xanderdir/bin/run_xander_skel.sh $Setenvfile "build" "$Genes"', file=fw_qsub)
        print('echo "*** build finished.."', file=fw_qsub)
        
        print('echo "*** find starting.."', file=fw_qsub)
        print('bash $Xanderdir/bin/run_xander_skel.sh $Setenvfile "find" "$Genes"', file=fw_qsub)
        print('echo "*** find finished.."', file=fw_qsub)
        
        print('echo "*** search starting.."', file=fw_qsub)
        print('bash $Xanderdir/bin/run_xander_skel.sh $Setenvfile "search" "$Genes"', file=fw_qsub)
        print('echo "*** search finished.."', file=fw_qsub)

if __name__ == '__main__':
    main()
