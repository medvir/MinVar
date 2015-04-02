#!/usr/local/bin/python3.4

import os
import os.path
import sys
import subprocess
import shutil
import warnings

from Bio import SeqIO

import callvar
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


def find_subtype(sampled_reads=200, remote=False):
    ''''''
    import pandas as pd

    loc_hit_file = 'loc_res.tsv'
    rem_hit_file = 'rem_res.tsv'
    # columns defined in standard blast output format 6
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    # sample reads with seqtk and convert to fasta
    s_cmline = 'seqtk sample high_quality.fastq {} | seqret -filter -out \
                sample_hq.fasta'.format(sampled_reads)
    sout = subprocess.check_output(s_cmline, shell=True,
                                   universal_newlines=True)

    # prepare blast to subject sequences
    blast_cmline_0 = 'blastn -task megablast -query sample_hq.fasta -outfmt 6'
    blast_cmline = blast_cmline_0 + ' -subject {} -out {}'.format(references_file,
                                                                  loc_hit_file)
    bout = subprocess.check_output(blast_cmline, shell=True,
                                   universal_newlines=True)

    # read blast results and assigns to best subjects
    loc_hits = pd.read_csv(loc_hit_file, names=cols, delimiter="\t")
    n_hits = loc_hits.shape[0]
    queries = len(set(loc_hits['qseqid']))
    print('Queries: {}\tHits: {}'.format(queries, n_hits))
    freqs = {k: 0.0 for k in set(loc_hits['sseqid'])}

    grouped = loc_hits.groupby(['qseqid'])
    for name, group in grouped:
        matching = group[group['pident'] == group['pident'].max()]['sseqid']
        # if query has more max matching, shares the weight
        for m in matching:
            freqs[m] += 1./len(matching)

    locfreqs = ((k, freqs[k]/queries) for k in sorted(freqs, key=freqs.get, reverse=True) if freqs[k])
    return locfreqs

    '''
    if best_hit.pident < 97.:
        warnings.warn('Low identity (%5.3f)' % best_hit.pident)
        alarm += 1
    if best_hit.qcovs < 95.:
        warnings.warn('Low coverage (%5.3f)' % best_hit.qcovs)
        alarm += 1
    if best_hit.length < 100:
        warnings.warn('Short match (%d)' % best_hit.length)
        alarm += 1
    '''
    if remote > 1:
        warnings.warn('Too many warnings: running remotely.')
        blast_cmline = 'blastn -task megablast -query sample_hq.fasta -db nr'
        blast_cmline += ' -remote -entrez_query \"hiv1[organism]\"'
        blast_cmline += ' -outfmt 6 -out remote_results.tsv'
        bout = subprocess.check_output(blast_cmline, shell=True,
                                       universal_newlines=True)
        cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        hits = pd.read_csv('remote_results.tsv', names=cols, delimiter="\t")
        best_hit = hits.iloc[0]
        return best_hit.sseqid
        
    return


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
    '''
    group2 = parser.add_argument_group('Run options',
                                       'Parameters that can (maybe should) be \
                                        changed according to the needs')

    group3 = parser.add_argument_group('More options',
                                       'Do you really want to change this?')
    '''
    group1.add_argument("-a", "--amplicon", default=def_amp, type=str, dest="a",
                        help="fasta file with the amplicon")

    group1.add_argument("-f", "--fastq", default="", type=str, dest="f",
                        help="input reads in fastq format")
    
    '''group2.add_argument("-s", "--winshifts", default=3, type=int, dest="s",
                        help="number of window shifts <%(default)s>")

    group2.add_argument("-r", "--region", default='', type=str, dest="r",
                        help="region in format 'chr:start-stop',\
                        eg 'ch3:1000-3000'")'''

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


def make_consensus(hxb2_ref):
    '''Take high quality reads, align with blast to HXB2, make consensus'''

    # take 400 fasta reads
    cml = 'seqtk seq -A high_quality.fastq | seqtk sample - 800 > hq_smp.fasta'
    subprocess.call(cml, shell=True)

    # blast them
    cml = 'blastn -task blastn -subject %s -query hq_smp.fasta -outfmt 5 \
            -out outblast.xml -word_size 7 -qcov_hsp_perc 80' % hxb2_ref
    subprocess.call(cml, shell=True)

    # convert to SAM -> BAM -> sort, index
    b2s_exe = os.path.join(dn, 'blast2sam.py')
    cml = b2s_exe
    cml += ' outblast.xml > outblast.sam'
    subprocess.call(cml, shell=True)

    cml = 'samtools view -F 16 -u outblast.sam | samtools sort - outblast_sorted'
    subprocess.call(cml, shell=True)

    cml = 'samtools index outblast_sorted.bam'
    subprocess.call(cml, shell=True)

    # mpileup and consensus with bcftools consensus (v1.2 required)
    cml = 'samtools mpileup -uf %s outblast_sorted.bam | bcftools call -mv -Oz -o calls.vcf.gz' % hxb2_ref
    subprocess.call(cml, shell=True)

    cml = 'tabix calls.vcf.gz'
    subprocess.call(cml, shell=True)

    cml = 'cat %s | bcftools consensus calls.vcf.gz > tmpcns.fasta' % hxb2_ref
    subprocess.call(cml, shell=True)

    # change sequence name to sample_cons, remove old file
    cml = 'sed -e \'s/HXB2_pol/sample_cons/\' tmpcns.fasta > cns.fasta'
    subprocess.call(cml, shell=True)    
    os.remove('tmpcns.fasta')


def main(read_file=None, amp_file=None):

    try:
        args = parse_com_line()
        amp_file = args.a
    except:
        amp_file = def_amp
    if not read_file:
        read_file = args.f
    assert os.path.exists(read_file), 'File %s not found' % read_file

    # locally defined filter_reads writes reads into high_quality.fastq
    filter_reads(read_file)

    # consensus is saved into cns.fasta
    make_consensus(amp_file)

    # find subtype: only used in the report
    try:
        matched_types = find_subtype()
    except:
        match_id = ''

    # align reads, call mutations and write the report
    al_reads = align_reads(ref='cns.fasta', reads='high_quality.fastq')
    mut_file = callvar.main(ref_file='cns.fasta', bamfile=al_reads)
    reportdrm.main(mut_file, matched_types)


if __name__ == "__main__":
    main()
