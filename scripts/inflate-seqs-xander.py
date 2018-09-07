#! /usr/bin/env python
# by gjr, 09012016

import screed
import sys
import os
import math

def main():
    if len(sys.argv) != 3:
        mes = ('*** python {} <infile.fa> <file.coverage>')
        print >> sys.stderr, mes.format(os.path.basename(sys.argv[0]))
        sys.exit(1)

    if sys.argv[1] == '-':
        input_fp = sys.stdin
    else:
        input_fp = open(sys.argv[1])

    d = {}
    for line in open(sys.argv[2]):
        if line.startswith('#'):
            continue
        lis = line.split('\t')
        name = lis[0]
        cov = int(math.ceil(float(lis[1])))
        d[name] = cov

    n = 0
    for n, record in enumerate(screed.fasta.fasta_iter(input_fp)):
        name = record.name
        seq = record.sequence
        if name == '#=GC_RF':
            print '>{}\n{}'.format(name, seq)
            continue
        if name not in d:
            continue
        cov = d[name]
        for i in range(cov):
            print '>{}__copy{}\n{}'.format(name, i, seq)

    n += 1

if __name__ == '__main__':
    main()
