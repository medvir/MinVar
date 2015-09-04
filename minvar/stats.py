#!/usr/local/bin/python3.4

import os
import sys
import subprocess

dn = os.path.dirname(__file__)
bed_file = os.path.join(dn, 'db/consensus_B.bed')

def coverage_stats(bam_file):
    '''Computes the coverage for the different genes by invoking
    samtools bedcov
    '''
    cml = 'samtools bedcov %s %s' % (bed_file, bam_file)
    proc = subprocess.Popen(cml, shell=True, stdout=subprocess.PIPE,
                            universal_newlines=True)
    with proc.stdout as handle:
        for l in handle:
            lsp = l.strip().split('\t')
            cov_per_base = float(lsp[-1]) / (int(lsp[2]) - int(lsp[1]) + 1)
            print(lsp, cov_per_base)

    cum_cov = {}
    cml = 'bedtools coverage -b %s -a %s -hist' % (bam_file, bed_file)
    proc = subprocess.Popen(cml, shell=True, stdout=subprocess.PIPE,
                            universal_newlines=True)
    with proc.stdout as handle:
        for l in handle:
            lsp = l.strip().split('\t')
            if len(lsp) < 8:
                continue
            region = lsp[3]
            cov_here = int(lsp[4])
            fract_at_cov = float(lsp[7])
            if cov_here > 100:
                cum_cov[region] = cum_cov.get(region, 0) + fract_at_cov
    for k, v in cum_cov.items():
        print(k, v)

if __name__ == "__main__" and __package__ is None:
    coverage_stats(sys.argv[1])
