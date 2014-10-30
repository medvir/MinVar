#!/usr/local/bin/python3.4

import os
import os.path
import sys
import subprocess
import shutil
import warnings

from Bio import SeqIO

import callvar
from OptimAssembly import optimassembly
from ReportDRM.src import reportdrm

import logging
import logging.handlers

# Make a global logging object.
mylog = logging.getLogger(__name__)
dn = os.path.dirname(__file__)

references_file = os.path.join(dn, 'db/HIV_consensus.fasta')
def_amp = os.path.join(dn, 'db/HXB2_pol_gene.fasta')
qual_thresh = 20
min_len = 75

def filter_reads(filename):
    '''Use seqtk and Biopython to trim and filter low quality reads'''

    # run seqtk trimfq to trim low quality ends
    mylog.info('Trimming reads with seqtk')
    if filename.endswith('gz'):
        mylog.info('Reads in gzip format')
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

    alarm = 0
    if best_hit.pident < 97.:
        warnings.warn('Low identity (%5.3f)' % best_hit.pident)
        alarm += 1
    if best_hit.qcovs < 95.:
        warnings.warn('Low coverage (%5.3f)' % best_hit.qcovs)
        alarm += 1
    if best_hit.length < 100:
        warnings.warn('Short match (%d)' % best_hit.length)
        alarm += 1
    if alarm > 1:
        sys.exit('Too many warnings: consider run remotely.')
    # TODO: write the alternative to run remotely
    return best_hit.name


def align_reads(ref=None, reads=None):
    '''Align all high quality reads to found reference with smalt,
    convert and sort with samtools'''

    cml = 'smalt index -k 7 -s 2 ref %s' % ref
    subprocess.call(cml, shell=True)
    cml = 'smalt map -n 12 -o hq_2_cons.sam -x -c 0.8 -y 0.8 ref %s' % reads
    subprocess.call(cml, shell=True)
    cml = 'samtools view -Su hq_2_cons.sam | samtools sort - hq_2_cons_sorted'
    subprocess.call(cml, shell=True)
    cml = 'samtools index hq_2_cons_sorted.bam'
    subprocess.call(cml, shell=True)

    os.remove('hq_2_cons.sam')
    return 'hq_2_cons_sorted.bam'


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

    group1.add_argument("-a", "--amplicon", default=def_amp, type=str, dest="a",
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


def main(read_file=None, amp_file=None):
    import subprocess

    if not read_file:
        args = parse_com_line()
        read_file = args.f
    if not amp_file:
        amp_file = def_amp

    assert os.path.exists(read_file), 'File %s not found' % read_file
    # locally defined filter_reads writes reads into high_quality.fastq
    filter_reads(read_file)
    ref = list(SeqIO.parse(amp_file, 'fasta'))[0]
    mylog.info('Calling optimassembly')
    optimassembly.main('high_quality.fastq', amp_file, len(ref), 2)
    
    # find subtype: only used in the report
    try:
        match_id = find_subtype()
        #al_ref = 'matching_reference.fasta'
        #cml = 'seqret %s:%s -outseq %s -auto' % (references_file, match_id, al_ref)
        #subprocess.call(cml, shell=True)
    except:
        match_id = ''

    hxb2_ref = amp_file
    # TODO hxb2ref should be inferred from the amplicon
    al_reads = align_reads(ref=hxb2_ref, reads='high_quality.fastq')
    mut_file = callvar.main(ref_file=amp_file, bamfile=al_reads)
    reportdrm.main(mut_file, match_id)


if __name__ == "__main__":
    main()