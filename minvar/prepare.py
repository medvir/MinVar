#!/usr/bin/env python3

import os
import os.path
import sys
import gzip
import shlex
import subprocess
import shutil
import warnings
import logging
from pkg_resources import resource_filename

from Bio import SeqIO

references_file = resource_filename(__name__, 'db/HIV_cons_db.nsq')
def_amp = resource_filename(__name__, 'db/consensus_B.fna')

qual_thresh = 20

try:
    CPUS = max(1, os.cpu_count())
except TypeError:
    CPUS = 1


def compute_min_len(filename):
    ''' Estimate distribution of sequence length and find a reasonable minimum
    length.
    '''
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    logging.info('Read input file to record lengths')
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as handle:
            reads_len = [len(s) for (n, s, q) in FastqGeneralIterator(handle)]
    else:
        with open(filename, 'r') as handle:
            reads_len = [len(s) for (n, s, q) in FastqGeneralIterator(handle)]

    reads_len = sorted(reads_len)
    n = len(reads_len)
    percentile_5 = reads_len[int(0.05 * n)]
    logging.info('5th percentile: %d' % percentile_5)
    return percentile_5 - 2


def filter_reads(filename, max_n, min_len=49):
    '''Use seqtk and Biopython to trim and filter low quality reads'''
    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    # run seqtk trimfq to trim low quality ends
    logging.info('Trimming reads with seqtk')
    r1 = 'seqtk trimfq %s | seqtk sample - %d' % (filename, max_n)
    oh = open('high_quality.fastq', 'w')
    proc = subprocess.Popen(r1, shell=True, stdout=subprocess.PIPE,
                            universal_newlines=True)
    with proc.stdout as handle:
        for name, s, quals in FastqGeneralIterator(handle):
            rl = len(s)
            if 'N' in s or rl < min_len:
                #float(sum(quals)) / rl < qual_thresh or
                continue
            else:
                print('@%s\n%s\n+\n%s' % (name, s, quals), file=oh)
    oh.close()
    return oh.name

def find_subtype(sampled_reads=1000, remote=False):
    ''''''
    import pandas as pd

    loc_hit_file = 'loc_res.tsv'
    rem_hit_file = 'rem_res.tsv'
    # columns defined in standard blast output format 6
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    # sample reads with seqtk and convert to fasta
    s_cmline = 'seqtk sample high_quality.fastq {} | seqtk seq -A - > \
                sample_hq.fasta'.format(sampled_reads)
    sout = subprocess.check_output(s_cmline, shell=True,
                                   universal_newlines=True)

    # prepare blast to subject sequences
    ref_stem = os.path.splitext(references_file)[0]
    blast_cmline = 'blastn -task megablast -query sample_hq.fasta -outfmt 6 -num_threads %d' % min(6, CPUS)
    blast_cmline += ' -db {} -out {}'.format(ref_stem, loc_hit_file)
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
            freqs[m] += 1. / len(matching)

    oh = open('subtype_evidence.csv', 'w')
    for k in sorted(freqs, key=freqs.get, reverse=True):
        if freqs[k]:
            print('%s,%5.4f' % (k, freqs[k]/queries), file=oh)
    oh.close()
    return
    #
    # if best_hit.pident < 97.:
    #     warnings.warn('Low identity (%5.3f)' % best_hit.pident)
    #     alarm += 1
    # if best_hit.qcovs < 95.:
    #     warnings.warn('Low coverage (%5.3f)' % best_hit.qcovs)
    #     alarm += 1
    # if best_hit.length < 100:
    #     warnings.warn('Short match (%d)' % best_hit.length)
    #     alarm += 1
    #
    # if remote:
    #     warnings.warn('Too many warnings: running remotely.')
    #     blast_cmline = 'blastn -task megablast -query sample_hq.fasta -db nr'
    #     blast_cmline += ' -remote -entrez_query \"hiv1[organism]\"'
    #     blast_cmline += ' -outfmt 6 -out remote_results.tsv'
    #     bout = subprocess.check_output(blast_cmline, shell=True,
    #                                    universal_newlines=True)
    #     cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    #             'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    #     hits = pd.read_csv('remote_results.tsv', names=cols, delimiter="\t")
    #     best_hit = hits.iloc[0]
    #     return best_hit.sseqid
    #
    # return


def align_reads(ref=None, reads=None, out_file=None, mapper='bwa'):
    '''Align all high quality reads to found reference with smalt,
    convert and sort with samtools'''

    if mapper == 'smalt':
        cml = 'smalt index -k 7 -s 2 cnsref %s' % shlex.quote(ref)
        subprocess.call(shlex.split(cml))
        cml = 'smalt map -n %d -o hq_2_cons.sam -x -c 0.8 -y 0.8 cnsref %s' % (min(12, CPUS), reads)
        subprocess.call(shlex.split(cml))
    elif mapper == 'bwa':
        cml = 'bwa index -p cnsref %s' % shlex.quote(ref)
        subprocess.call(shlex.split(cml))
        cml = 'bwa mem -t %d -O 12 cnsref %s > hq_2_cons.sam' % (min(12, CPUS), reads)
        subprocess.call(cml, shell=True)
    elif mapper == 'novo':
        cml = 'novoindex cnsref.ndx %s' % ref
        subprocess.call(cml, shell=True)
        cml = 'novoalign -d cnsref.ndx -f %s -F STDFQ -o SAM > hq_2_cons.sam' % reads
        subprocess.call(cml, shell=True)

    cml = 'samtools view -Su hq_2_cons.sam | samtools sort -T /tmp -@ %d -o %s -' % (min(4, CPUS), out_file)
    subprocess.call(cml, shell=True)
    cml = 'samtools index %s' % out_file
    subprocess.call(cml, shell=True)

    os.remove('hq_2_cons.sam')

    return out_file


def phase_variants(reffile, varfile):
    '''Parsing variants in varfile (vcf) and applying them to reffile
    '''

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_dna

    wobble = {
        'AG' : 'R',
        'CT' : 'Y',
        'AC' : 'M',
        'GT' : 'K',
        'CG' : 'S',
        'AT' : 'W'
    }

    logging.info('Phasing file %s with variants in %s' % (reffile, varfile))
    refseq = list(SeqIO.parse(reffile, 'fasta'))[0]
    reflist = list(str(refseq.seq))
    reflist.insert(0, '')

    with open(varfile) as h:
        for l in h:
            if l.startswith('#'):
                continue
            lsp = l.split('\t')
            pos = int(lsp[1])
            ref, alt = lsp[3:5]
            infos = dict(a.split('=') for a in lsp[7].split(';'))
            varlen = len(ref)
            assert str(refseq.seq)[pos - 1:pos - 1 + varlen] == ref, '%s instead of %s at pos %d' % (str(refseq.seq)[pos - 1:pos - 1 + varlen], ref, pos)
            # exclude multiallelic positions
            nalleles = len(alt.split(','))
            assert nalleles == 1
            # reference and alternate frequency
            freqs = 1. - float(infos['AF']), float(infos['AF'])

            if freqs[0] < 0.5:
                report = alt
            else:
                report = ref
            for i, r in enumerate(report):
                reflist[pos + i] = r

    finalrefseq = ''.join(reflist)
    fs = SeqRecord(Seq(finalrefseq), id='sample_cons_Pol', description='lofreq')
    SeqIO.write(fs, 'tmp_cns.fasta', 'fasta')
    return 'tmp_cns.fasta'


def make_consensus(ref_file, reads_file, out_file, sampled_reads=4000,
                   mapper='bwa', cons_caller='own'):
    '''Take reads, align to reference, return consensus file'''
    import glob
    import time
    import shlex

    pf = time.perf_counter()
    ranseed = int(str(pf).split('.')[1][-5:])

    # sample reads
    cml = 'seqtk sample -s %d %s %d > hq_smp.fastq' % (ranseed, reads_file, sampled_reads)
    subprocess.call(cml, shell=True)
    logging.info('reads sampled')

    if mapper == 'blast':
        logging.info('take fasta reads')
        cml = 'seqtk seq -A hq_smp.fastq > hq_smp.fasta'
        subprocess.call(cml, shell=True)

        logging.info('blast them')
        cml = 'blastn -task blastn -subject %s -query hq_smp.fasta -outfmt 5\
                -out outblast.xml -word_size 7 -qcov_hsp_perc 80' % ref_file
        subprocess.call(cml, shell=True)

        logging.info('convert to SAM -> BAM')
        cml = 'blast2sam.py outblast.xml > refcon.sam'
        subprocess.call(cml, shell=True)

        # reverse reads are not yet properly treated, so -F 16
        cml = 'samtools view -F 16 -u refcon.sam | samtools sort -T /tmp -o refcon_sorted.bam'
        subprocess.call(cml, shell=True)

    elif mapper == 'bwa':
        logging.info('mapping reads with bwa')
        cml = 'bwa index -p tmpref %s' % shlex.quote(ref_file)
        subprocess.call(shlex.split(cml))

        cml = shlex.split('bwa mem -t %d -O 12 tmpref hq_smp.fastq' % min(CPUS, 12))
        with open('refcon.sam', 'w') as f:
            subprocess.call(cml, stdout=f)

        logging.debug('remove index files')
        for biw in glob.glob('tmpref.*'):
            os.remove(biw)
        logging.debug('SAM -> BAM -> SORT')
        cml1 = shlex.split('samtools view -u refcon.sam')
        view = subprocess.Popen(cml1, stdout=subprocess.PIPE)
        cml2 = shlex.split('samtools sort -T /tmp -@ %d -o refcon_sorted.bam' % min(6, CPUS))
        sort = subprocess.check_output(cml2, stdin=view.stdout)
        view.wait()

    elif mapper == 'novo':
        logging.info('mapping with novoalign')
        cml = 'novoindex tmpref.ndx %s' % ref_file
        subprocess.call(cml, shell=True)
        cml = 'novoalign -c %d -d tmpref.ndx -f hq_smp.fastq -F STDFQ -o SAM > refcon.sam' % min(12, CPUS)
        subprocess.call(cml, shell=True)
        os.remove('tmpref.ndx')

    #os.remove('refcon.sam')
    cml = shlex.split('samtools index refcon_sorted.bam')
    subprocess.call(cml)

    if cons_caller == 'bcftools':
        logging.info('mpileup and consensus with bcftools consensus (v1.2 required)')
        cml = 'samtools mpileup -uf %s refcon_sorted.bam | bcftools call -mv -Ov -o calls.vcf' % ref_file
        subprocess.call(cml, shell=True)

        var_pres = False
        with open('calls.vcf') as ch:
            for l in ch:
                if not l.startswith('#'):
                    logging.ionfo('variants are present')
                    var_pres = True
                    break

        if var_pres:
            cml = 'bgzip -f calls.vcf'
            subprocess.call(cml, shell=True)
            cml = 'tabix calls.vcf.gz'
            subprocess.call(cml, shell=True)
            cml = 'cat %s | bcftools consensus calls.vcf.gz > tmp_cns.fasta' % ref_file
            subprocess.call(cml, shell=True)
            os.remove('calls.vcf.gz.tbi')
        else:
            shutil.copyfile(ref_file, 'tmp_cns.fasta')

    elif cons_caller == 'own':
        logging.info('applying variants to consensus')
        # we need to copy the ref_file locally because lofreq parallel does
        # not quote filenames and fails with sapces/parentheses
        shutil.copy(ref_file, '.')
        local_ref = os.path.split(ref_file)[1]

        cml = shlex.split('samtools faidx %s' % local_ref)
        subprocess.call(cml)

        cml = 'lofreq call-parallel --pp-threads %d -f %s refcon_sorted.bam -o calls.vcf' % (min(CPUS, 6), local_ref)
        subprocess.call(shlex.split(cml))

        phase_variants(ref_file, 'calls.vcf')
        cml = 'bgzip -f calls.vcf'
        subprocess.call(cml, shell=True)

    os.rename('tmp_cns.fasta', out_file)
    return out_file


def iterate_consensus(reads_file, ref_file=def_amp):
    ''' Call make_consensus until convergence or for a maximum number of
    iterations
    '''

    logging.info('First consensus iteration')
    iteration = 1
    # consensus from first round is saved into cns_1.fasta
    cns_file_1 = make_consensus(ref_file, reads_file,
                                out_file = 'cns_%d.fasta' % iteration,
                                sampled_reads=50000, mapper='bwa')
    try:
        os.rename('calls.vcf.gz', 'calls_%d.vcf.gz' % iteration)
    except FileNotFoundError:
        loggin.info('No variants found in making consensus')
    dh = compute_dist('cns_%d.fasta' % iteration, def_amp)
    logging.info('distance is %6.4f' % dh)

    if dh < 5:
        logging.info('Converged at first step')
        return cns_file_1

    while iteration <= 10:
        logging.info('iteration %d' % (iteration + 1))
        # use cns_1.fasta for further consensus rounds
        new_cons = make_consensus('cns_%d.fasta' % iteration, reads_file,
                                  out_file='cns_%d.fasta' % (iteration + 1),
                                  sampled_reads=50000, mapper='bwa')
        try:
            os.rename('calls.vcf.gz', 'calls_%d.vcf.gz' % (iteration + 1))
        except FileNotFoundError:
            loggin.info('No variants found in making consensus')
        assert new_cons == 'cns_%d.fasta' % (iteration + 1)
        iteration += 1
        dh = compute_dist('cns_%d.fasta' % iteration, new_cons)
        logging.info('distance is %6.4f' % dh)
        if dh < 5:
            logging.info('converged!')
            break

    return new_cons


def compute_dist(file1, file2):
    cml = 'blastn -query %s -subject %s -out pair.tsv -outfmt "6 qseqid sseqid pident"' % (shlex.quote(file1), shlex.quote(file2))

    subprocess.call(cml, shell=True)
    with open('pair.tsv') as h:
        for l in h:
            ident = l.split('\t')[-1]
    return 100 - float(ident)


def main(read_file=None, max_n_reads=200000):
    assert os.path.exists(read_file), 'File %s not found' % read_file

    min_len = compute_min_len(read_file)

    # locally defined filter_reads writes reads into high_quality.fastq
    filtered_file = filter_reads(read_file, max_n_reads, min_len)

    find_subtype()

    cns_file = iterate_consensus(filtered_file)

    # change sequence name to sample_cons, remove old file
    cml = 'sed -e \'s/CONSENSUS_B/sample_cons_Pol/\' cns_2.fasta > cns_final.fasta'
    subprocess.call(cml, shell=True)
    os.remove('cns_2.fasta')
    cml = 'samtools faidx cns_final.fasta'
    subprocess.call(cml, shell=True)


    prepared_file = align_reads(ref='cns_final.fasta', reads=filtered_file, out_file='hq_2_cns_final.bam')

    return 'cns_final.fasta', prepared_file

if __name__ == "__main__":
    main()
