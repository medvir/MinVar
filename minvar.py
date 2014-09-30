#!/usr/local/bin/python3.4

import os
import os.path
import sys
import subprocess
import shutil

from Bio import SeqIO

import logging
import logging.handlers

# Make a global logging object.
mylog = logging.getLogger(__name__)

dn = os.path.dirname(__file__)

references_file = os.path.expanduser('HIV_consensus.fasta')

qual_thresh = 20
min_len = 75
from OptimAssembly import optimassembly

def filter_reads(filename):
    '''Use seqtk and Biopython to trim and filter low quality reads'''

    # run seqtk trimfq to trim low quality ends
    mylog.info('Trimming reads with seqtk')
    if filename.endswith('gz'):
        optlog.info('Reads in gzip format')
        r1 = 'gunzip -c %s | seqtk trimfq - ' % filename
    else:
        r1 = 'seqtk trimfq %s' % filename

    oh = open('high_quality.fastq', 'w')
    proc = subprocess.Popen(r1, shell=True, stdout=subprocess.PIPE,
                            universal_newlines=True)
    with proc.stdout as handle:
        for s in SeqIO.parse(handle, 'fastq'):
            quals = s.letter_annotations['phred_quality']
            rl = len(s)
            if 'N' in s or \
                float(sum(quals)) / rl < qual_thresh or \
                rl < min_len:
                continue
            else:
                SeqIO.write(s, oh, 'fastq')
    oh.close()


def find_subtype():
    ''''''
    import pandas as pd

    blast_out_file = 'blast_cons_ref.tsv'
    format_specs = ['qseqid', 'sseqid', 'score', 'pident', 'qcovs',
                       'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                       'sstart', 'send']
    blast_out = open(blast_out_file, 'w')
    blast_out.write('\t'.join(format_specs) + '\n')
    blast_cmline = 'blastn -task megablast -query consensus.fasta'
    blast_cmline += ' -subject %s' % references_file
    blast_cmline += ' -outfmt \'6 %s\'' % ' '.join(format_specs)
    bout = subprocess.check_output(blast_cmline, shell=True,
                                   universal_newlines=True)
    blast_out.write(bout)
    blast_out.close()
    hits = pd.read_csv(blast_out_file, index_col='sseqid', delimiter="\t")
    hits.sort(['score', 'pident', 'length', 'qcovs'], ascending=False,
              inplace=True)
    best_hit = hits.iloc[0]
    if best_hit.pident < 98. or best_hit.qcovs < 95. or \
        best_hit.length < 100:
        sys.exit('Good hit not found: run remotely!')
    # TODO: write the alternative to run remotely
    return best_hit.name


def align_reads(ref=None, reads=None):
    '''Align all high quality reads to found reference with smalt,
    convert and sort with samtools'''

    cml = 'smalt index -k 7 -s 2 consref %s' % ref
    subprocess.call(cml, shell=True)
    cml = 'smalt map -n 12 -o hq_2_cons.sam -x -c 0.8 -y 0.8 consref %s' % reads
    subprocess.call(cml, shell=True)
    cml = 'samtools view -Su hq_2_cons.sam | samtools sort - hq_2_cons_sorted'
    subprocess.call(cml, shell=True)
    cml = 'samtools index hq_2_cons_sorted.bam'
    subprocess.call(cml, shell=True)

    os.remove('hq_2_cons.sam')


def parse_com_line():

    import argparse

    # parse command line
    parser = argparse.ArgumentParser()

    # First define all option groups
    group1 = parser.add_argument_group('Input files', 'Required input')

    group2 = parser.add_argument_group('Run options',
                                       'Parameters that can (maybe should) be \
                                        changed according to the needs')

    group3 = parser.add_argument_group('More options',
                                       'Do you really want to change this?')

    group1.add_argument("-a", "--amplicon", default='', type=str, dest="a",
                        help="fasta file with the amplicon")

    group1.add_argument("-f", "--fastq", default="", type=str, dest="f",
                        help="input reads in fastq format")

    group2.add_argument("-s", "--winshifts", default=3, type=int, dest="s",
                        help="number of window shifts <%(default)s>")

    group3.add_argument("-i", "--sigma", default=0.01, type=float, dest="i",
                        help="value of sigma to use when calling\
                        SNVs <%(default)s>")

    group2.add_argument("-r", "--region", default='', type=str, dest="r",
                        help="region in format 'chr:start-stop',\
                        eg 'ch3:1000-3000'")

    # set logging level
    mylog.setLevel(logging.DEBUG)
    # This handler writes everything to a file.
    LOG_FILENAME = './minvar.log'
    h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w',
                                             maxBytes=100000, backupCount=5)
    f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s\
                          %(lineno)d %(message)s")
    h.setFormatter(f)
    mylog.addHandler(h)
    mylog.info(' '.join(sys.argv))

    return parser.parse_args()


if __name__ == "__main__":

    import subprocess

    args = parse_com_line()

    assert os.path.exists(args.f), 'File %s not found' % args.f
    # locally defined filter_reads writes reads into high_quality.fastq
    #filter_reads(args.f)
    ref = list(SeqIO.parse(args.a, 'fasta'))[0]
#    optimassembly.main('high_quality.fastq', args.a, len(ref), 2)
    # find subtype and use it for alignment 
    match_id = find_subtype()
    al_ref = 'matching_reference.fasta'
    cml = 'seqret %s:%s -outseq %s -auto' % (references_file, match_id, al_ref)
    subprocess.call(cml, shell=True)
    align_reads(ref=al_ref, reads='high_quality.fastq')
    