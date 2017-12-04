#!/usr/bin/env python3
"""A few functions to obtain statistics on the sample."""
import sys
import shlex
import subprocess


def coverage_stats_per_base(bam_file, bed_file):
    """Compute the coverage for the different genes by invoking samtools bedcov."""
    cml = shlex.split(
        'samtools bedcov %s %s' % (shlex.quote(bed_file), bam_file))
    proc = subprocess.Popen(
        cml, stdout=subprocess.PIPE, universal_newlines=True)
    with proc.stdout as handle:
        for l in handle:
            lsp = l.strip().split('\t')
            cov_per_base = float(lsp[-1]) / (int(lsp[2]) - int(lsp[1]) + 1)
            print(lsp, cov_per_base)


def gene_coverage(bam_file, bed_file, threshold=100):
    """Compute fraction of genes covered by more than threshold reads."""
    cum_cov = {}
    cml = shlex.split('bedtools coverage -b %s -a %s -hist' %
                      (bam_file, shlex.quote(bed_file)))
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
            if cov_here > threshold:
                cum_cov[region] = cum_cov.get(region, 0) + fract_at_cov
    return cum_cov


def genome_coverage(bam_file, threshold=5):
    """Invoke bedtools genomecov to detect zero coverage regions."""
    # bedtools genomecov -ibam hq_2_cns_final_recal.bam -bga | grep -w 0$
    cml1 = shlex.split('bedtools genomecov -ibam %s' % bam_file)
    bed = subprocess.run(cml1, stdout=subprocess.PIPE,
                         universal_newlines=True)
    cov_fract = 0.0
    cov_pos = 0
    for l in bed.stdout.splitlines():
        if not l.startswith('genome'):
            continue
        lsp = l.split()
        depth, positions, fract = int(lsp[1]), int(lsp[2]), float(lsp[4])
        if depth >= threshold:
            cov_fract += fract
            cov_pos += positions
    return cov_fract, cov_pos


def genome_longest_covered(bam_file, threshold=20):
    """Invoke bedtools to find the longest/most covered contiguous region."""
    import os
    import pandas as pd
    # The command used is:
    # bedtools genomecov -ibam refcon_sorted.bam -bg | awk '$4 > threshold' | bedtools cluster -i -
    # giving a file with lines whose meaning is
    # chr          start    stop  coverage   cluster
    # e.g.
    # CONSENSUS_B     56      59      101     1

    # WATCH OUT: this seems to be 0-based

    # convert coverage into bedgraph format, only keep regions above threshold
    cml = shlex.split('bedtools genomecov -ibam %s -bg' % bam_file)
    proc = subprocess.Popen(cml, stdout=subprocess.PIPE, universal_newlines=True)
    with proc.stdout as reading, open('covfile.txt', 'w+') as writing:
        for l in reading:
            lsp = l.strip().split('\t')
            cov = int(lsp[3])
            if cov >= threshold:
                writing.write(l)

    # cluster adiacent entries in covfile
    cml = shlex.split('bedtools cluster -i covfile.txt')
    with open('clusters.txt', 'w+') as oh:
        subprocess.call(cml, stdout=oh, universal_newlines=True)
    # read them in a dataframe, compute the lenght of each region and return the longest
    df = pd.read_table('clusters.txt', header=None, names=['chr', 'start', 'stop', 'coverage', 'cluster'])
    clusters = df.groupby('cluster').agg({'start':'min', 'stop':'max'}).reset_index()
    clusters['length'] = clusters['stop'] - clusters['start']
    best = clusters.loc[clusters['length'].idxmax()]
    os.remove('covfile.txt')
    os.remove('clusters.txt')
    return best.to_dict()


if __name__ == "__main__" and __package__ is None:
    # coverage_stats_per_base(sys.argv[1], sys.argv[2])
    glg = genome_longest_covered(sys.argv[1], 400)
    print(glg)
