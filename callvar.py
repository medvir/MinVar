#!/usr/local/bin/python3.4

import sys
import subprocess
import os
import warnings
import pandas as pd
import numpy as np
from pprint import pprint
from ReportDRM.src import Alignment

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

RAW_DEPTH_THRESHOLD = 100
MIN_FRACTION = 0.025
MAPPING_QUALITY_THRESHOLD = 50

# amminoacid sequences from files in db directory
dn_dir = os.path.dirname(__file__)
db_dir = os.path.abspath(os.path.join(dn_dir, 'ReportDRM/db'))
prot = \
    list(SeqIO.parse(os.path.join(db_dir, 'protease.faa'), 'fasta'))[0]
RT = list(SeqIO.parse(os.path.join(db_dir, 'RT.faa'), 'fasta'))[0]
integrase = \
    list(SeqIO.parse(os.path.join(db_dir, 'integrase.faa'), 'fasta'))[0]

# 64 codons + '---'
translation_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*', '---': '-'}

info_fields = {}


def find_frame(ref):
    '''Returns frame [1, 2, 3] based on length of longest ORF.
    It uses gapped translation with the table defined locally, since Biopython
    does not support it (yet?)'''

    max_len = -1
    refseq = str(list(SeqIO.parse(ref, 'fasta'))[0].seq)
    refseq = refseq.replace('-', '').upper()

    if len(refseq) < 100:
        warnings.warn('Estimation of the frame not reliable: ' +
                      'sequence too short')
    aa_seq = {}
    for f in [1, 2, 3]:
        nts = refseq[f - 1:]
        codons = [nts[i:i + 3] for i in range(0, len(nts), 3)]
        aa_seq[f] = ''.join([translation_table.get(c, '*') for c in codons])

        # get maximal length of ORFs
        len1 = max([len(a) for a in aa_seq[f].split('*')])
        if len1 > max_len:
            best_frame = f
            max_len = len1

    return best_frame, refseq[best_frame - 1:], aa_seq[best_frame]


def filter_variant(info_field):
    '''Parses the info field in vcf file and returns false if one filter is
    triggered'''

    infos = dict(a.split('=') for a in info_field.split(';'))

    # Alternate allele observations
    freqs = [float(f) for f in infos['AO'].split(',')]

    # if filters are not met, put allele observation to zero

    # alternate observations on forward strand
    for i, info in enumerate(infos['SAF'].split(',')):
        if int(info) < 5:
            freqs[i] = 0.0
            warnings.warn('SAF filter triggered')
    # alternate observations on reverse strand
    for i, info in enumerate(infos['SAR'].split(',')):
        if int(info) < 5:
            freqs[i] = 0.0
            warnings.warn('SAR filter triggered')
    # reads placed left supporting alternate
    for i, info in enumerate(infos['RPL'].split(',')):
        if int(info) < 5:
            freqs[i] = 0.0
            warnings.warn('RPL filter triggered')
    # reads placed right supporting alternate
    for i, info in enumerate(infos['RPR'].split(',')):
        if int(info) < 5:
            freqs[i] = 0.0
            warnings.warn('RPR filter triggered')
    # strand balance probablity for alternate obs
    for i, info in enumerate(infos['SAP'].split(',')):
        if float(info) < 20:
            freqs[i] = 0.0
            warnings.warn('SAP filter triggered')

    # insert observations of the reference allele at the beginning and filter
    freqs.insert(0, float(infos['RO']))
    if float(infos['SRP']) < 10 or int(infos['SRR']) < 5 or int(infos['SRF']) < 5:
        warnings.warn('Putting ref frequency to zero')
        print('Freq 2 0')
        freqs[0] = 0.0

    # lost observations should be few
    if abs(1. - sum(freqs) / int(infos['DP'])) > 0.1:
        pprint('Lost observations: DP=%s, RO+AO=%d' % (infos['DP'], sum(freqs)))

    # returns normalised frequencies of reference, allele1 [allele2, ...]
    return [f / int(infos['DP']) for f in freqs]


def parsevar(frame, nt_ref, aa_ref, vcf_file):
    '''Parses mutations in vcf file and returns them wrt to consensus B.
    In order to achieve this several steps are needed. Mutations are called by
    freebayes wrt to consensus of reads. This function writes sample consensus
    sequences modified by the variants and compare them to consensus B. It then
    returns the mutations wrt consensus B avoiding the duplication of those
    already found on the sample consensus.
    '''
    import re

    # nt_ref and aa_ref are already framed
    Dn = 0
    Ds = 0
    allvars = []

    prog = re.compile('##INFO=<ID=(.*),Number.*,Description="(.*)">')
    all_muts = []
    vc_lines = []
    # pass on the file, saving info fields too
    with open(vcf_file) as h:
        for l in h:
            if l.startswith('#'):
                res = prog.match(l)
                try:
                    info_fields[res.group(1)] = res.group(2)
                    continue
                except:
                    continue
            vc_lines.append(l)

    vcf_mutations = pd.DataFrame(columns=['wt', 'pos', 'mut', 'freq', 'gene'])
    # vcf_mutations = vcf_mutations.set_index(['gene', 'pos', 'mut'])
    # now parses the variants
    for l in vc_lines:
        lsp = l.split('\t')
        # firs of all, filter on quality
        if float(lsp[5]) < 100:
            continue

        # pos is wrt sample consensus, check that the reference matches
        pos = int(lsp[1])
        ref_nt = lsp[3]
        assert lsp[3] == nt_ref[pos - frame:pos - frame + len(ref_nt)], 'ref shift: %s' % lsp

        # pos - frame is 0-based position on the framed nt reference
        # consequently, (pos - frame) / 3 is on aa reference
        codon_pos = int((pos - frame) / 3)

        # triplet will be on positions (codon_pos * 3, (codon_pos + 1) * 3)
        # check that translation matches
        wt_cod = nt_ref[codon_pos * 3: (codon_pos + 1) * 3]
        assert translation_table[wt_cod] == aa_ref[codon_pos]

        # here several filters are implemented, freqs is populated with the
        # normalised frequencies of [ref, allele1 [, allele2, allele3...] ]
        try:
            freqs = filter_variant(lsp[7])
        except:
            print(lsp)
            print("Unexpected error:", sys.exc_info()[0])
            sys.exit()
        if not freqs:
            continue

        tmp_freq = freqs[0]

        # save alt alleles; loop needed when there is more than one
        alt = lsp[4]  # alt can contain more than one allele
        for i, alt_nt in enumerate(alt.split(',')):
            mut_nt = nt_ref[:pos - frame] + alt_nt
            mut_nt += nt_ref[pos - frame + len(alt_nt):]
            mut_seq = Seq(mut_nt, generic_dna)
            mut_aa = str(mut_seq.translate())
            # mut_aa will already contain mutations found in sample consensus
            # plus the (translated) mutation alt_nt found in the vcf file.
            # Only the latter must be called wrt consensus B
            len_aa_var = len(alt_nt) / 3
            for mreg in [prot, RT, integrase]:
                pm_here = parse_mutations(str(mreg.seq), mut_aa, mreg.id,
                                          freqs[i + 1],
                                          codon_pos, codon_pos + len_aa_var)
                if pm_here:
                    vcf_mutations = vcf_mutations.append(pm_here)
            # save to fasta file
            sr = SeqRecord(mut_seq.translate(),
                           id='%s-%s-%s;freq=%f' % (lsp[3], lsp[1], lsp[4],
                                                    freqs[i + 1]),
                           description='')
            # dict_here = {'gene': gene_name, 'wt': z[0],
            #             'pos': i, 'mut': z[1], 'freq': freq_here}
            tmp_freq += freqs[i + 1]
            all_muts.append(sr)
        assert tmp_freq <= 1.0, ','.join(freqs) + '<--->' + str(tmp_freq)

    SeqIO.write(all_muts, 'mutants.fasta', 'fasta')
    # print('Dn: %d\tDs: %d' % (Dn, Ds))
    return vcf_mutations


def parse_mutations(ref_seq, input_seq, gene_name, freq_here,
                    start=None, stop=None):
    '''Takes amino acids sequence and aligns against a reference to compute
    mutations, returns a panda data frame.
    start and stop can be used to give the restrict to a region on input_seq'''

    dict_list = []
    Alignment.needle_align('asis:%s' % ref_seq,
                           'asis:%s' % input_seq, 'h.tmp')
    alhr = Alignment.alignfile2dict(['h.tmp'])
    os.remove('h.tmp')
    alih = alhr['asis']['asis']
    alih.summary()
    if 3 * alih.mismatches > alih.ident:
        warnings.warn('Too many mismatches')
    wt, mut = alih.seq_a, alih.seq_b

    i = 0
    end_wt, end_mut = len(wt.rstrip('-')), len(mut.rstrip('-'))
    end_pos = min(end_wt, end_mut)
    mut_pos = 0
    for z in list(zip(wt, mut))[:end_pos]:
        i += z[0] != '-'
        mut_pos += z[1] != '-'
        if z[0] != z[1] and i:
            mutation = '%s%d%s' % (z[0], i, z[1])
            dict_here = {'gene': gene_name, 'wt': z[0],
                         'pos': i, 'mut': z[1], 'freq': freq_here}
            if (not start and not stop) or (start <= mut_pos and mut_pos <= stop):
                dict_list.append(dict_here)
    return dict_list


def main(ref_file='cns.fasta', bamfile='hq_2_cons_sorted.bam', parallel=True, n_regions=8):
    '''What does a main do?'''

    # parse mutations already found in the consensus wrt to
    # consensus B prot, RT and integrase
    frame, nt_framed, aa_framed = find_frame(ref_file)
    cns_mutations = pd.DataFrame()
    for mreg in [prot, RT, integrase]:
        cns_mutations = cns_mutations.append(parse_mutations(str(mreg.seq),
                                             aa_framed, mreg.id, 1.0))
    cns_mutations.rename(columns={'freq': 'freq_cns'}, inplace=True)
    cns_mutations.drop('wt', axis=1, inplace=True)
    cns_mutations.to_csv('cns_mutations.csv', index=False)
    bamstem = '.'.join(bamfile.split('.')[:-1])

    # call minority variants with freebayes
    if parallel:
        # first compute regions
        cml = 'samtools faidx %s' % ref_file
        subprocess.call(cml, shell=True)
        cml = '/usr/local/freebayes/bamtools/bin/bamtools coverage -in %s > xyz' % bamfile
        subprocess.call(cml, shell=True)
        ctr_exe = os.path.join(dn_dir, 'coverage_to_regions.py')
        cml = 'cat xyz | %s %s.fai %d > regions.bed' % (ctr_exe, ref_file, n_regions)
        subprocess.call(cml, shell=True)
        os.remove('xyz')
        # now run freebayes-parallel
        # need to change PATH first
        os.environ["PATH"] += os.pathsep + '/usr/local/freebayes/vcflib/bin/'
        cml = '/usr/local/freebayes/scripts/freebayes-parallel regions.bed %d' % n_regions
    else:
        cml = 'freebayes'

    cml += ' --min-alternate-count 15 --min-coverage %s' % RAW_DEPTH_THRESHOLD
    cml += ' --min-alternate-fraction %f --ploidy 1 --pooled-continuous' % MIN_FRACTION
    cml += ' --haplotype-length 30'
    cml += ' --fasta-reference %s %s > %s_fb.vcf' % (ref_file, bamfile, bamstem)
    if os.path.exists('%s_fb.vcf' % bamstem):
        print('vcf file exists, reading from it', file=sys.stderr)
    else:
        print('freebayes will compute variants', file=sys.stderr)
        subprocess.call(cml, shell=True)

    # parse mutations in vcf file
    vcf_mutations = parsevar(frame, nt_framed, aa_framed, '%s_fb.vcf' % bamstem)

    # sometimes the reference is put to 0 frequency in parsevar, might be
    # useful later
    zero_freqs = vcf_mutations[vcf_mutations['freq'] == 0.0]

    vcf_mutations = vcf_mutations[vcf_mutations['freq'] > 0.0]
    # synonymous mutations can be merged
    vcf_mutations = vcf_mutations.groupby(['gene', 'pos', 'mut']).sum().reset_index()
    vcf_mutations.to_csv('vcf_mutations.csv', index=False)
    # manipulate mutations
    out = pd.merge(cns_mutations, vcf_mutations, on=['gene', 'pos', 'mut'], how='outer')
    minfreq = out.filter(like='freq').min(axis=1, skipna=True)
    out.drop(['freq', 'freq_cns'], axis=1, inplace=True)
    out['freq'] = minfreq

    cols = ['gene', 'pos', 'mut', 'freq']
    out.sort(cols, inplace=True)
    out = out[cols]
    out['pos'] = out['pos'].astype(int)
    out['freq'] = np.round(out['freq'], 4)
    out.to_csv('mutations.csv', sep=',', index=False)
    return 'mutants.fasta'

if __name__ == '__main__':
    main()
