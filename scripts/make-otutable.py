import os
import subprocess

def parse_mcclust(f, target_cutoff):
    """
    Convert parse mcclust.clust to dictionary of cluster# and seqname list

    Parameters:
    -----------
    f : str
        clustering result (.clust file) from mcclust
    target_cutoff: str
        distance cutoff used for OTU (e.g., 0.03)

    Returns:
    --------
    a dictionary with cluster# as key and sequence name set as value

    """
    
    fp = open(f)
    target_cutoff = float(target_cutoff)
    assert 0 <= target_cutoff <= 1
    triger = False

    d = {}
    temp_cutoff = None
    temp_total = None
    temp_str = None
    for line in fp:
        if 'File' in line:
            continue
        if 'Sequences:' in line:
            continue
        line = line.strip()
        if 'distance cutoff:' in line:
            _cutoff = line.split(':',1)[1].strip()
            continue
        if 'Total Clusters:' in line:
            total = line.split(':', 1)[1].strip()
            total = int(total)
            triger = True
            continue

        if triger:
            if not line:
                _cutoff = float(_cutoff)

                if _cutoff == target_cutoff:
                    return d
                elif _cutoff > target_cutoff:
                    return temp_d

                temp_cutoff = _cutoff
                temp_total = total
                temp_d = d
                d = {}
                triger = False
                continue

            assert len(line.split('\t')) == 4, 'parsing wrong ..'
            cluNum, s, num, names = line.split('\t')
            cluNum = int(cluNum)
            name_set = set(names.split())
            for name in name_set:
                d[name] = cluNum


def main():
    SAMPLES = snakemake.config['SAMPLES']
    JAR_DIR = os.path.abspath(snakemake.config['JAR_DIR'])
    MAX_JVM_HEAP = snakemake.config['MAX_JVM_HEAP'] # memory for java program

    aligned_files = [os.path.abspath(f) for f in snakemake.input.alignedprot_list]
    coverage_files = [os.path.abspath(f) for f in snakemake.input.covfile_list]
    outdir=os.path.abspath(os.path.dirname(snakemake.output[0]))

    dist=snakemake.config['DISTANCE']  # range 0 to 0.5 

    #derep = "cd {} && java -Xmx{} -jar {}/Clustering.jar derep -o derep.fa -m '#=GC_RF' ids samples {} || {{ echo 'derep failed' ;  exit 1; }}; "
    #dmatrix = "cd {} && java -Xmx{} -jar {}/Clustering.jar dmatrix  -c 0.5 -I derep.fa -i ids -l 50 -o dmatrix.bin || {{ echo 'dmatrix failed' ;  exit 1; }}; "
    #cluster = "cd {} && java -Xmx{} -jar {}/Clustering.jar cluster -d dmatrix.bin -i ids -s samples -o complete.clust || {{ echo 'cluster failed' ;  exit 1; }}; "
    #otutable = "cat {} > merged_covfile && java -Xmx{} -jar {}/Clustering.jar cluster_to_Rformat complete.clust {} {} {} merged_covfile; "

    #cmd = (
    #    "mkdir -p {}; ".format(outdir)
    #    + derep.format(outdir, MAX_JVM_HEAP, JAR_DIR, ' '.join(aligned_files))
    #    + dmatrix.format(outdir, MAX_JVM_HEAP, JAR_DIR)
    #    + cluster.format(outdir, MAX_JVM_HEAP, JAR_DIR)
    #    + "cd {} && rm dmatrix.bin nonoverlapping.bin; ".format(outdir)
    #)

    derep = "java -Xmx{} -jar {}/Clustering.jar derep -o derep.fa -m '#=GC_RF' ids samples {} || {{ echo 'derep failed' ;  exit 1; }}; "
    dmatrix = "java -Xmx{} -jar {}/Clustering.jar dmatrix  -c 0.5 -I derep.fa -i ids -l 50 -o dmatrix.bin || {{ echo 'dmatrix failed' ;  exit 1; }}; "
    cluster = "java -Xmx{} -jar {}/Clustering.jar cluster -d dmatrix.bin -i ids -s samples -o complete.clust || {{ echo 'cluster failed' ;  exit 1; }}; "

    cmd = (
        "mkdir -p {}; cd {}; ".format(outdir, outdir)
        + derep.format(MAX_JVM_HEAP, JAR_DIR, ' '.join(aligned_files))
        + dmatrix.format(MAX_JVM_HEAP, JAR_DIR)
        + cluster.format(MAX_JVM_HEAP, JAR_DIR)
        + "rm dmatrix.bin nonoverlapping.bin; ".format(outdir)
    )

    subprocess.run(cmd, shell=True, check=True)

    # make otutable from mcclust output
    target_cutoff = dist
    clust_listfile = '{}/complete.clust'.format(outdir)
    outfile = '{}/otutable.tsv'.format(outdir)
    cov_files = coverage_files

    d = parse_mcclust(clust_listfile,target_cutoff)
    with open(outfile, 'w') as fw:
        otu_list = ['OTU{}'.format(otu) for otu in sorted(set(d.values()))]
        otu_num = len(otu_list)
        print('{}\t{}'.format('Sample', '\t'.join(otu_list)), file=fw)
        for cov_f in cov_files:
            d_otu_cov = {}
            with open(cov_f) as fp:
                for line in fp:
                    if line.startswith('#'):
                        continue
                    line = line.rstrip()
                    #seqid  mean_cov        median_cov      total_pos       covered_pos     covered_ratio
                    _lis = line.split()
                    name = _lis[0]
                    cov = float(_lis[1])

                    otu = d[name]
                    d_otu_cov[otu] = d_otu_cov.get(otu, 0) + cov

            for i in range(1, otu_num+1):
                d_otu_cov[i] = d_otu_cov.get(i, 0)

            # sort by otu
            items = sorted(d_otu_cov.items())
            assert len(d_otu_cov) == otu_num
            list = ['{:.1f}'.format(abu) for otu, abu in items]
            tag = os.path.basename(cov_f).split('.')[0].split('_')[0]
            print('{}\t{}'.format(tag, '\t'.join(list)), file=fw)


if __name__ == '__main__':
    main()
