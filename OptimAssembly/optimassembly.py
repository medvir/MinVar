#!/usr/bin/env python

import sys
import random
import subprocess

import logging
import logging.handlers

from Bio import SeqIO
from collections import Counter

# Make a global logging object.
optlog = logging.getLogger(__name__)

replicates = 5
qual_thresh = 30
exp_cov = 100
min_len = 75


def print_fasta(reads, filename):
    import textwrap

    oh = open(filename, 'w')
    for i, rep_r in enumerate(reads):
        oh.write('>read_%d\n' % i)
        tx = textwrap.fill(rep_r, 80)
        oh.write(tx + '\n')
    oh.close()


def filter_reads(filename):
    '''Use seqtk and Biopython to trim and filter low quality reads'''

    # run seqtk trimfq to trim low quality ends
    optlog.info('Trimming reads with seqtk')
    if filename.endswith('gz'):
        optlog.info('Reads in gzip format')
        r1 = 'gunzip -c %s | seqtk trimfq - ' % filename
    else:
        r1 = 'seqtk trimfq %s' % filename

    c = 1
    sequence_bag = set([])
    oh = open('selected.fastq', 'w')
    cnt = Counter()
    proc = subprocess.Popen(r1, shell=True, stdout=subprocess.PIPE,
                            universal_newlines=True)
    optlog.info('Filtering')
    with proc.stdout as handle:
        for s in SeqIO.parse(handle, 'fastq'):
            quals = s.letter_annotations['phred_quality']
            rl = len(s)
            seq = str(s.seq)
            if 'N' in seq or \
                float(sum(quals)) / rl < qual_thresh or \
                rl < min_len:
                continue
            if seq in sequence_bag:
                SeqIO.write(s, oh, 'fastq')
                c += 1
                cnt.update([seq])
                if c % 50000 == 0:
                    print >> sys.stderr, 'Written %10d reads' % c
            else:
                sequence_bag.add(seq)
    oh.close()
    return cnt


def compute_msa(ref):
    '''Takes contigs in velvet output directory, checks the orientation and
    computes MSA
    '''
    save_contigs = []
    for rep in range(1, replicates + 1):
        contigs_file = 'out_%d/contigs.fa' % rep
        contigs = [s for s in SeqIO.parse(contigs_file, 'fasta')]
        optlog.info('In replicate %d found %d contigs' % (rep, len(contigs)))
        for cont in contigs:
            # forward
            exe_string = 'needle %s asis:%s -auto -out stdout | grep Score | cut -d \" \" -f 3' % \
                (ref, cont.seq)
            score = subprocess.check_output(exe_string, shell=True)
            score_forward = float(score)
            # reverse
            exe_string = 'needle %s asis:%s -sreverse2 -auto -out stdout | grep Score | cut -d \" \" -f 3' % \
                (ref, cont.seq)
            score = subprocess.check_output(exe_string, shell=True)
            score_reverse = float(score)

            if score_forward >= score_reverse:
                save_contigs.append(cont)
            else:
                save_contigs.append(cont.reverse_complement())

    allc = SeqIO.write(save_contigs, 'all_contigs.fasta', 'fasta')
    optlog.info('All %d contigs written to all_contigs.fasta' % allc)

    optlog.info('Trying multiple alignment with muscle')
    exe_string = 'muscle -maxiters 1 -diags -in all_contigs.fasta -out msa.fasta &> log_muscle.log'
    msa = subprocess.call(exe_string, shell=True)
    if msa > 1:
        optlog.error('Error %d returned trying muscle' % msa)
        return
    else:
        optlog.info('muscle run, see output in msa.fasta')

    optlog.info('Trying to get consensus')
    exe_string = 'cons msa.fasta -outseq consensus.fasta -plurality 0.0001 -auto \
    -name consensus_contigs'
    msa = subprocess.call(exe_string, shell=True)
    if msa > 1:
        optlog.error('Error %d returned trying cons' % msa)
        return
    else:
        optlog.info('cons run, see output in consensus.fasta')

    optlog.info('Trying alignment with the reference')
    try:
        exe_string = 'needle %s consensus.fasta  -auto -aformat3 sam -stdout | grep -v ^@ | cut -f 6' % ref
        cigar = subprocess.check_output(exe_string, shell=True)
        optlog.info('Alignment run, cigar string is %s' % cigar)
        print >> sys.stderr, 'Reference vs. consensus cigar: %s' % cigar
    except CalledProcessError:
        optlog.error('Could not run the alignment')


def run_velvet(fastafile, dir_n, contig_cut):
    exe_string = 'velveth out_%d 21 -short -fasta %s' % (dir_n, fastafile)
    optlog.info('running' + exe_string)
    subprocess.call(exe_string, shell=True)
    exe_string = 'velvetg out_%d -unused_reads yes -cov_cutoff auto \
                  -exp_cov auto -min_contig_lgth %d' % (dir_n, contig_cut)
    optlog.info('running' + exe_string)
    subprocess.call(exe_string, shell=True)


def main(filename, reference, length):
    '''
    '''

    # set logging level
    optlog.setLevel(logging.DEBUG)
    # This handler writes everything to a file.
    LOG_FILENAME = './optimas.log'
    hl = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w',
                                              maxBytes=100000, backupCount=5)
    f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s\
                          %(lineno)d %(message)s")
    hl.setFormatter(f)
    optlog.addHandler(hl)
    optlog.info(' '.join(sys.argv))

    # if average quality is low, or Ns are present, or read is short discard
    cnt = filter_reads(filename)

    optlog.info('There are %d unique reads' % len(cnt))

    non_singletons = [p[0] for p in cnt.items() if p[1] != 1]
    optlog.info(' %d are not singletons' % len(non_singletons))
    # integer division used in mean length
    mean_len = sum((len(st) for st in non_singletons))
    mean_len /= len(non_singletons)
    optlog.info('%d is their mean length' % mean_len)
    out_n_reads = length * exp_cov / mean_len
    optlog.info('%d needed for the desired coverage %d' %
                (out_n_reads, exp_cov))

    for r in range(replicates):
        try:
            optlog.info('Taking random sample %d' % r)
            rep_reads = random.sample(non_singletons, out_n_reads)
            fname = 'replicate_%d.fasta' % (r + 1)
            print_fasta(rep_reads, fname)
            run_velvet(fname, r + 1, mean_len)
        except ValueError:
            optlog.info('Taking all reads')
            fname = 'non_singleton_reads.fasta'
            print_fasta(non_singletons, fname)
            run_velvet(fname, r, mean_len)
            break

    # only if reference is given and replicates are present
    if reference and r:
        compute_msa(reference)


if __name__ == "__main__":

    import argparse
    # parse command line
    parser = argparse.ArgumentParser(description='Optimise de novo assembly\
                                     for short, viral genomes')
    parser.add_argument('-f', '--fastq', dest='fastq', type=str,
                        default='',
                        help='input file in fastq format <%(default)s>')
    parser.add_argument('-r', '--reference', dest='reference', type=str,
                        help='closest known genome reference')
    parser.add_argument('-l', '--length', dest='exp_length', type=int,
                        default=10000, help='expected length <%(default)s>')

    args = parser.parse_args()
    main(filename=args.fastq, reference=args.reference, length=args.exp_length)
