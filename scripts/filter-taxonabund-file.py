import sys
import os

def main():
    lis_infile = snakemake.input.taxonabunfile_list
    lis_tag = snakemake.config['SAMPLES']
    outfile = snakemake.output[0]

    lis_pair = zip(lis_infile, lis_tag)
    with open(outfile, 'w') as fw:
        print('Sample\tDomain\tPhylum\tGenus\tMatchName\tAbun\tFrac', file=fw)
        for infile, tag in lis_pair:
            triger = False
            for line in open(infile):
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('Lineage\tMatchName'):
                    triger = True
                    continue
                if triger:
                    str_taxa, str_rest  = line.split('\t', 1)
                    lis_taxa = str_taxa.split(';')

                    if len(lis_taxa) < 3:
                        lis_taxa = lis_taxa + ['Other']*(3-len(lis_taxa))

                    domain = lis_taxa[0]
                    if 'environmental' in domain.lower():
                        domain = 'Other'

                    phylum = lis_taxa[1]
                    if phylum == 'Proteobacteria':
                        phylum = lis_taxa[2]
                    if 'environmental' in phylum.lower():
                        phylum = 'Other'

                    if len(lis_taxa) <= 4:
                        genus = 'Other'
                    else:
                        genus = lis_taxa[-1]

                    if 'environmental' in genus.lower():
                        genus = 'Other'

                    new_line = '{}\t{}\t{}\t{}\t{}'.format(tag, domain, 
                                                 phylum, genus,
                                                 str_rest)

                    print(new_line, file=fw)

if __name__ == '__main__':
    main()
