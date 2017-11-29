#!/usr/bin/env python3
"""Prepare reads: filter, find organism and subtype, make consensus."""
import gzip
import logging
import os
import os.path
import shlex
import shutil
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pkg_resources import resource_filename

from .common import hcv_map, org_dict
from .stats import (genome_coverage, genome_longest_covered)

dn_dir = os.path.dirname(os.path.abspath(__file__))
HCV_references = resource_filename(__name__, 'db/HCV/subtype_references.fasta')
HIV_references = resource_filename(__name__, 'db/HIV/subtype_references.fasta')
HCV_recomb_references = \
    resource_filename(__name__, 'db/HCV/recomb_references.fasta')
HIV_recomb_references = \
    resource_filename(__name__, 'db/HIV/recomb_references.fasta')
blast2sam_exe = 'blast2sam'  # resource_filename(__name__, '../scripts/blast2sam.py')
HIV_amp = resource_filename(__name__, 'db/HIV/consensus_B.fna')
# HCV_amp =

qual_thresh = 20

try:
    CPUS = max(1, os.cpu_count())
except TypeError:
    CPUS = 1


def disambiguate(dna_string):
    """Removes ambiguous bases from a sequence"""
    from random import choice
    d2a = {'AG': 'R', 'CT': 'Y', 'AC': 'M', 'GT': 'K', 'CG': 'S', 'AT': 'W'}
    wobbles = {v: k for k, v in d2a.items()}
    out_seq = [s if s not in wobbles.keys() else choice(wobbles[s]) for s in dna_string]
    return ''.join(out_seq)
# def pad_consensus(denovo_seq, organism, subtype):
#     """denovo consensus will most likely not span the whole sequenced region.
#
#     The rest of the program expects sequences of a specific length, so we pad
#     the denovo sequence with the missing region.
#     """
#     if organism == 'HIV':  # use the mutation neutral consensus_B
#         ref = list(SeqIO.parse(HIV_amp, 'fasta'))[0]
#     elif organism == 'HCV':  # use the best subtype
#         if subtype == 'RF1_2k/1b':
#             ref = list(SeqIO.parse(HCV_recomb_references, 'fasta'))[0]
#             assert ref.id.startswith('AY587845')
#         else:
#             an = acc_numbers[subtype][0]
#             ref = [s for s in SeqIO.parse(HCV_references, 'fasta')
#                    if s.id.startswith(an)][0]
#     al_file = 'padalign.fasta'
#     needle_align('asis:%s' % str(ref.seq),
#                  'asis:%s' % str(denovo_seq.seq), al_file, go=20.0, ge=4.0)
#
#     al_ref, al_den = [str(s.seq) for s in
#                       SeqIO.parse(al_file, 'fasta')]
#     #os.remove(al_file)
#     # first non gap positions
#     den_start = len(al_den) - len(al_den.lstrip('-'))
#     den_stop = len(al_den.rstrip('-'))
#     ref_start = len(al_ref) - len(al_ref.lstrip('-'))
#     ref_stop = len(al_ref.rstrip('-'))
#     # where the shortest sequence ends
#     min_stop = min(den_stop, ref_stop)
#     # if denovo doesn't cover the beginning, take it from the reference
#     if den_start > ref_start:
#         padded = al_ref[:den_start] + al_den[den_start:min_stop]
#     else:
#         padded = al_den[ref_start:min_stop]
#     # if denovo doesn't cover the end, take it from the reference
#     if ref_stop > min_stop:
#         padded += al_ref[min_stop:ref_stop]
#
#     return padded


def compute_min_len(filename):
    """Estimate distribution of sequence length and find a reasonable minimum length."""
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
    logging.info('5th percentile: %d', percentile_5)
    return percentile_5 - 2


def filter_reads(filename, max_n, min_len=49):
    """Use seqtk and Biopython to trim and filter low quality reads."""
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
                # float(sum(quals)) / rl < qual_thresh or
                continue
            else:
                print('@%s\n%s\n+\n%s' % (name, s, quals), file=oh)
    oh.close()
    return oh.name


def find_subtype(reads_file, sampled_reads=1000, recomb=False):
    """Blast a subset of reads against references to infer organism and support for group/subtype."""
    import fileinput
    import pandas as pd

    if recomb:
        ref_files = [HIV_recomb_references, HCV_recomb_references]
    else:
        ref_files = [HIV_references, HCV_references]

    sub_ref = 'temp_ref.fasta'
    with open(sub_ref, 'w') as fout, fileinput.input(ref_files) as fin:
        for line in fin:
            fout.write(line)

    loc_hit_file = 'loc_res.tsv'
    # columns defined in standard blast output format 6
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    # sample reads with seqtk and convert to fasta, equivalent to:
    # seqtk sample reads_file n_reads | seqtk seq -A - > sample_hq.fasta
    cml1 = shlex.split('seqtk sample %s %d' % (reads_file, sampled_reads))
    sample = subprocess.Popen(cml1, stdout=subprocess.PIPE)
    cml2 = shlex.split('seqtk seq -A -')
    with open('sample_hq.fasta', 'w') as oh:
        to_fasta = subprocess.Popen(cml2, stdin=sample.stdout, stdout=oh)
        sample.stdout.close()
        output = to_fasta.communicate()[0]
    del output
    oh.close()

    # prepare blast to subject sequences
    blast_cmline = 'blastn -task megablast -query sample_hq.fasta -outfmt 6'
    # -num_threads %d' % min(6, CPUS)
    blast_cmline += ' -subject {} -out {}'.format(shlex.quote(sub_ref),
                                                  loc_hit_file)
    blast = subprocess.Popen(shlex.split(blast_cmline), universal_newlines=True)
    blast.wait()
    os.remove(sub_ref)

    # read blast results and assigns to best subjects
    loc_hits = pd.read_csv(loc_hit_file, names=cols, delimiter="\t")
    n_hits = loc_hits.shape[0]
    queries = len(set(loc_hits['qseqid']))
    print('Queries: {}\tHits: {}'.format(queries, n_hits))
    freqs = {k: 0.0 for k in set(loc_hits['sseqid'])}
    support = {'HIV': 0.0, 'HCV': 0.0}
    grouped = loc_hits.groupby(['qseqid'])
    for name, group in grouped:
        del name
        # take the best matching hit, i.e. the one with max percent identity
        matching = group[group['pident'] == group['pident'].max()]['sseqid']
        # if query has more max matching, shares the weight
        for m in matching:
            freqs[m] += 1. / (queries * len(matching))
            support[org_dict[m.split('.')[0]]] += 1. / len(matching)

    max_freq = max(freqs.values())
    organism = max(support, key=support.get)
    freq2 = {}
    for k, v in freqs.items():
        if organism == 'HCV':
            # replace accession numbers with genotypes
            gt = hcv_map[k.split('.')[0]]
        else:
            # HIV references already have subtypes in the name
            gt = k
        freq2[gt] = freq2.get(gt, 0.0) + freqs[k]
        # save sequence id of best hit
        if v == max_freq:
            max_acc = k

    return organism, freq2, max_acc

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
    """Align all high quality reads to found reference with smalt, convert and sort with samtools."""
    if mapper == 'smalt':
        cml = 'smalt index -k 7 -s 2 cnsref %s' % shlex.quote(ref)
        subprocess.call(shlex.split(cml))
        cml = 'smalt map -n %d -o hq_2_cons.sam -x -c 0.8 -y 0.8 cnsref %s' %\
            (min(12, CPUS), reads)
        subprocess.call(shlex.split(cml))
    elif mapper == 'bwa':
        cml = 'bwa index -p cnsref %s' % shlex.quote(ref)
        subprocess.call(shlex.split(cml))
        cml = 'bwa mem -t %d -O 12 cnsref %s > hq_2_cons.sam' %\
            (min(12, CPUS), reads)
        subprocess.call(cml, shell=True)
    elif mapper == 'novo':
        cml = 'novoindex cnsref.ndx %s' % ref
        subprocess.call(cml, shell=True)
        cml = 'novoalign -d cnsref.ndx -f %s -F STDFQ -o SAM > hq_2_cons.sam' \
            % reads
        subprocess.call(cml, shell=True)

    cml = \
        'samtools view -Su hq_2_cons.sam | samtools sort -T /tmp -@ %d -o %s -'\
        % (min(4, CPUS), out_file)
    subprocess.call(cml, shell=True)
    cml = 'samtools index %s' % out_file
    subprocess.call(cml, shell=True)
    os.remove('hq_2_cons.sam')

    return out_file


def phase_variants(reffile, varfile):
    """Parse variants in varfile (vcf) and applying them to reffile."""
    # from Bio.Alphabet import generic_dna
    logging.info('Phasing file %s with variants in %s', reffile, varfile)
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
            assert str(refseq.seq)[pos - 1:pos - 1 + varlen] == ref, \
                '%s instead of %s at pos %d' \
                % (str(refseq.seq)[pos - 1:pos - 1 + varlen], ref, pos)
            # exclude multiallelic positions
            nalleles = len(alt.split(','))
            assert nalleles == 1
            # reference and alternate frequency
            freqs = 1. - float(infos['AF']), float(infos['AF'])
            # do not put gaps in the consensus
            if alt == '-':
                report = ref
            elif ref == '-':
                report = alt
            else:
                if freqs[0] < 0.5:
                    report = alt
                else:
                    report = ref
            for i, r in enumerate(report):
                reflist[pos + i] = r

    finalrefseq = ''.join(reflist)
    fs = Seq(finalrefseq)
    return fs


def make_consensus(ref_file, reads_file, out_file, sampled_reads=10000,
                   mapper='bwa', cons_caller='own'):
    """Take reads, align to reference, return consensus file."""
    import glob
    from time import perf_counter

    pf = perf_counter()
    ranseed = int(str(pf).split('.')[1][-5:])

    # sample reads
    cml = shlex.split('seqtk sample -s %d %s %d'
                      % (ranseed, reads_file, sampled_reads))
    with open('hq_smp.fastq', 'w') as oh:
        subprocess.call(cml, stdout=oh)
    logging.info('reads sampled')

    if mapper == 'blast':
        logging.info('take fasta reads')
        cml = shlex.split('seqtk seq -A hq_smp.fastq')
        with open('hq_smp.fasta', 'w') as oh:
            subprocess.call(cml, stdout=oh)

        logging.info('blast them')
        cml = shlex.split('blastn -task blastn -subject %s -query hq_smp.fasta\
                          -outfmt 5 -out outblast.xml -word_size 7\
                          -qcov_hsp_perc 80' % ref_file)
        subprocess.call(cml)

        logging.info('convert to SAM -> BAM')
        cml = shlex.split('%s outblast.xml' % shlex.quote(blast2sam_exe))
        with open('refcon.sam', 'w') as oh:
            subprocess.call(cml, stdout=oh)

        # reverse reads are not yet properly treated, so -F 16
        cmls = 'samtools view -F 16 -u refcon.sam | samtools sort -T /tmp -o refcon_sorted.bam'
        subprocess.call(cmls, shell=True)

    elif mapper == 'bwa':
        logging.info('mapping reads with bwa')
        cml = shlex.split('bwa index -p tmpref %s' % shlex.quote(ref_file))
        subprocess.call(cml)

        cml = shlex.split('bwa mem -t %d -O 12 tmpref hq_smp.fastq' %
                          min(CPUS, 12))
        with open('refcon.sam', 'w') as f:
            subprocess.call(cml, stdout=f)

        logging.debug('remove index files')
        for biw in glob.glob('tmpref.*'):
            os.remove(biw)
        logging.debug('SAM -> BAM -> SORT')
        cml1 = shlex.split('samtools view -u refcon.sam')
        view = subprocess.Popen(cml1, stdout=subprocess.PIPE)
        cml2 = shlex.split('samtools sort -T /tmp -@ %d -o refcon_sorted.bam' %
                           min(6, CPUS))
        sort = subprocess.Popen(cml2, stdin=view.stdout,
                                stdout=subprocess.DEVNULL)
        view.stdout.close()
        output = sort.communicate()[0]

    os.remove('refcon.sam')
    cml = shlex.split('samtools index refcon_sorted.bam')
    subprocess.call(cml)
    covered_fract, covered_pos = genome_coverage('refcon_sorted.bam')
    logging.info('%d positions covered', covered_pos)

    if cons_caller == 'bcftools':
        logging.info(
            'mpileup and consensus with bcftools consensus (v1.2 required)')
        cml = 'samtools mpileup -uf %s refcon_sorted.bam | bcftools call -mv -Ov -o calls.vcf' % ref_file
        subprocess.call(cml, shell=True)

        var_pres = False
        with open('calls.vcf') as ch:
            for l in ch:
                if not l.startswith('#'):
                    logging.info('variants are present')
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
        # not quote filenames and fails with spaces/parentheses
        try:
            shutil.copy(ref_file, '.')
        except shutil.SameFileError:
            pass
        local_ref = os.path.split(ref_file)[1]

        cml = shlex.split('samtools faidx %s' % local_ref)
        subprocess.call(cml)

        cml = 'lofreq call-parallel --pp-threads %d -f %s refcon_sorted.bam -o calls.vcf' % (min(CPUS, 6), local_ref)
        subprocess.call(shlex.split(cml))

        phased_seq = phase_variants(ref_file, 'calls.vcf')
        cml = shlex.split('bgzip -f calls.vcf')
        subprocess.call(cml)

    phased_rec = SeqRecord(phased_seq, id='sample_consensus', description='')
    SeqIO.write(phased_rec, out_file, 'fasta')
    # del output
    return out_file, covered_fract


def iterate_consensus(reads_file, ref_file):
    """Call make_consensus until convergence or for a maximum number of iterations."""
    logging.info('First consensus iteration')
    iteration = 1
    # consensus from first round is saved into cns_1.fasta
    cns_file_1, new_cov = make_consensus(
        ref_file, reads_file, out_file='cns_%d.fasta' % iteration,
        sampled_reads=5000, mapper='blast')
    try:
        os.rename('calls.vcf.gz', 'calls_%d.vcf.gz' % iteration)
    except FileNotFoundError:
        logging.info('No variants found in making consensus')
    dh = compute_dist('cns_%d.fasta' % iteration, ref_file)
    logging.info('distance = %6.4f percent', dh)
    logging.info('covered = %6.4f ', new_cov)

    if dh < 2:
        logging.info('Converged at first step')
        return cns_file_1

    while iteration <= 10:
        logging.info('iteration %d', iteration + 1)
        # use cns_1.fasta for further consensus rounds

        new_cons, new_cov = make_consensus(
            'cns_%d.fasta' % iteration, reads_file,
            out_file='cns_%d.fasta' % (iteration + 1),
            sampled_reads=10000, mapper='bwa')
        try:
            os.rename('calls.vcf.gz', 'calls_%d.vcf.gz' % (iteration + 1))
        except FileNotFoundError:
            logging.info('No variants found in making consensus')
        assert new_cons == 'cns_%d.fasta' % (iteration + 1)
        iteration += 1
        dh = compute_dist('cns_%d.fasta' % iteration, new_cons)
        logging.info('distance is %6.4f', dh)
        logging.info('covered is %6.4f', new_cov)
        if dh < 2:
            logging.info('converged!')
            break

    return new_cons


def compute_dist(file1, file2):
    """Compute the distance in percent between two sequences (given in two files)."""
    cml = shlex.split(
        'blastn -gapopen 20 -gapextend 2 -query %s -subject %s -out pair.tsv \
         -outfmt "6 qseqid sseqid pident gaps"' % (shlex.quote(file1), shlex.quote(file2))
        )
    subprocess.call(cml)
    with open('pair.tsv') as h:
        for l in h:
            ident = l.split('\t')[2]
    os.remove('pair.tsv')
    return 100 - float(ident)


def main(read_file=None, max_n_reads=200000):
    """What the main does."""
    assert os.path.exists(read_file), 'File %s not found' % read_file

    min_len = compute_min_len(read_file)

    # locally defined filter_reads writes reads into high_quality.fastq
    filtered_file = filter_reads(read_file, max_n_reads, min_len)

    # infer organism and subtype
    organism, support_freqs, acc = find_subtype(filtered_file)
    logging.info('%s sequences detected', organism)
    max_support = max(support_freqs.values())

    # check if there is better support for recombinant (HCV only)
    max_support_rec = 0.0
    if organism == 'HCV':
        # if support is good, don't even try recombinants
        if max_support < 0.85:
            logging.info('Low support in HCV: try recombinant')
            organism, rec_support_freqs, rec_acc = find_subtype(filtered_file, recomb=True)
            max_support_rec = max(rec_support_freqs.values())
        # max_support_rec is 0.0, unless explicitely computed
        if max_support_rec > max_support:
            logging.info('Using recombinant')
            sub_file = HCV_recomb_references
            frequencies = rec_support_freqs
            s_id = rec_acc
        elif max_support_rec <= max_support:
            logging.info('Using non recombinant')
            sub_file = HCV_references
            frequencies = support_freqs
            s_id = acc
    elif organism == 'HIV':
        sub_file = HIV_references
        frequencies = support_freqs
        s_id = acc

    sorted_freqs = sorted(frequencies, key=frequencies.get, reverse=True)
    best_subtype = sorted_freqs[0]

    with open('subtype_evidence.csv', 'w') as oh:
        for k in sorted_freqs:
            if frequencies[k]:
                print('%s,%5.4f' % (k, frequencies[k]), file=oh)
    logging.info('Looking for best reference in file %s', sub_file)
    ref_dict = SeqIO.to_dict(SeqIO.parse(sub_file, 'fasta'))
    ref_rec = SeqRecord(ref_dict[s_id].seq, id=best_subtype.split('.')[0], description='')
    SeqIO.write([ref_rec], 'subtype_ref.fasta', 'fasta')
    cns_file = iterate_consensus(filtered_file, 'subtype_ref.fasta')
    # extract the longest region covered by at least 100 reads and save that
    longest_covered = genome_longest_covered('refcon_sorted.bam')
    all_ref = list(SeqIO.parse(cns_file, 'fasta'))[0]
    covered_dna = str(all_ref.seq[longest_covered['start'] - 1:longest_covered['stop']])
    all_ref.seq = Seq(disambiguate(covered_dna))
    SeqIO.write(all_ref, 'cns_final.fasta', 'fasta')
    #os.rename(cns_file, 'cns_final.fasta')
    #denovo_seq = list(SeqIO.parse('denovo_consensus.fasta', 'fasta'))[0]
    #padded = pad_consensus(denovo_seq, organism, best_subtype)
    #SeqIO.write(SeqRecord(Seq(padded), id='padded_consensus', description=''), 'cns_final.fasta', 'fasta')
    cml = shlex.split('samtools faidx cns_final.fasta')
    subprocess.call(cml)

    prepared_file = align_reads(ref='cns_final.fasta', reads=filtered_file, out_file='hq_2_cns_final.bam')

    return 'cns_final.fasta', prepared_file, organism


if __name__ == "__main__":
    import sys
    main(sys.argv[1])
