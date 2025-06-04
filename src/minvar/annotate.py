#!/usr/bin/env python3
"""
``annotate`` - Interprets the vcf and translates it into amminoacid mutations.

It is called by ``cli.py`` with the following parameters:

    :param vcf_file: VCF file with mutations
    :param ref_file: reference in FASTA format used to call mutations
    :param bam_file: BAM file that was used to call mutations
    :param organism: string defining what virus is being analyzed
    :type organism: HCV or HIV


``annotate.py`` is the most complex module in MinVar. By reading mutation and reference sequence
files, it builds a representation of all bases found at every position with relative frequency.
Further, it translates this representation into ammino acids and, aligning to organism references,
annotates them with gene name and position.

The main output is the file ``final.csv`` reporting all mutations with respect to the organism
consensus (consensus B for HIV, H77 for HCV). Example::

    gene,wt,pos,mut,freq
    E2,V,155,I,1.0000
    E2,T,175,S,1.0000
    E2,V,187,I,0.0270
    E2,V,191,A,0.9212
    E2,L,197,H,1.0000

Other two files produced are:

- ``cns_ambiguous.fasta`` nucleotide sequence, with wobble bases if frequency > 15%,
- ``cns_max_freq.fasta`` nucleotide sequence, only the nucleotide at max frequency is reported.

"""
import logging
import os
import shlex
import subprocess
import sys
import warnings

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

dn_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if __name__ == '__main__':
    if __package__ is None:
        os.sys.path.insert(1, dn_dir)
        mod = __import__('minvar')
        sys.modules["minvar"] = mod
        from common import h77_map, consensus_B_map, d2a, coverage_threshold, AMBIGUITY_threshold
        from Alignment import needle_align, alignfile2dict
else:
    from .common import h77_map, consensus_B_map, d2a, coverage_threshold, AMBIGUITY_threshold
    from .Alignment import needle_align, alignfile2dict


RAW_DEPTH_THRESHOLD = 50
VARIANT_QUALITY_THRESHOLD = 80
MAPPING_QUALITY_THRESHOLD = 20
HAPLO_FREQ_THRESHOLD = 0.015

# nt and amminoacid sequences from files in db directory
B_pol_nt_seq = list(
    SeqIO.parse(os.path.join(dn_dir, 'minvar/db/HIV/consensus_B.fna'), 'fasta')
)[0].seq
B_pol_aa_seq = list(
    SeqIO.parse(os.path.join(dn_dir, 'minvar/db/HIV/consensus_B.faa'), 'fasta')
)[0].seq
utr_5p_seq = list(
    SeqIO.parse(os.path.join(dn_dir, 'minvar/db/HCV/hcv_utr5.fasta'), 'fasta')
)[0].seq
h77_aa_seq = list(
    SeqIO.parse(os.path.join(dn_dir, 'minvar/db/HCV/H77_polyprotein.faa'),
                'fasta'))[0].seq
h77_nt_seq = list(
    SeqIO.parse(os.path.join(dn_dir, 'minvar/db/HCV/H77_cds.fna'),
                'fasta'))[0].seq

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
    """Return frame [1, 2, 3] based on length of longest ORF.

    It uses gapped translation with the table defined locally, since
    Biopython does not support it (yet?).
    """
    max_len = -1
    # refseq = str(ref)
    refseq = str(ref).replace('-', '').upper()
    assert not refseq.startswith('-')

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
    # to_trim = (len(refseq) - best_frame + 1) % 3
    return best_frame, aa_seq[best_frame].rstrip('*')


def parsevar(vcf_file, ref_seq):
    """Parse mutations in vcf file.

    This returns mutations in a pandas DataFrame.
    """
    import re
    # pass on the file, saving info fields too
    prog = re.compile('##INFO=<ID=(.*),Number.*,Description="(.*)">')
    vc_lines = []
    with open(vcf_file) as h:
        for l in h:
            if l.startswith('#'):
                res = prog.match(l)
                try:
                    info_fields[res.group(1)] = res.group(2)
                    continue
                except AttributeError:
                    continue
            vc_lines.append(l)

    # ref_seq = list(SeqIO.parse(ref_file, 'fasta'))[0]
    # frame, nt_framed, aa_framed = find_frame(ref_seq)
    ref_nt = str(ref_seq)
    vcf_mutations = []
    # vcf_mutations = vcf_mutations.set_index(['gene', 'pos', 'mut'])
    # now parses the variants
    for vc_line in vc_lines:
        lsp = vc_line.strip().split('\t')
        if float(lsp[5]) < VARIANT_QUALITY_THRESHOLD:
            logging.warning('Quality filter on position %s: qual=%s', lsp[1], lsp[5])
            continue

        pos = int(lsp[1])
        wt_nt, alt = lsp[3:5]
        try:
            assert ref_nt[pos - 1] == wt_nt[0], 'ref shift: %s not %s at pos %s' % (ref_nt[pos - 1], wt_nt, pos)
        except IndexError:
            print(ref_nt[pos - 10:pos + 10])
            print('problem at pos %s, wt_nt:%s' % (pos, wt_nt))
            sys.exit()
        infos = {}
        for a in lsp[7].split(';'):
            asp = a.split('=')
            try:
                infos[asp[0]] = asp[1]
            except IndexError:
                infos[asp[0]] = None

        # loop because multiple variants can be found on the same vcf line
        for alt_nt, alt_f in zip(alt.split(','), infos['AF'].split(',')):
            if len(alt_nt) > len(wt_nt):  # insertion
                logging.info('insertion found at pos %d', pos)
                # pm_here = {'wt': wt_nt, 'pos': pos, 'mut': alt_nt, 'freq': float(alt_f)}
                # vcf_mutations = vcf_mutations.append(pm_here, ignore_index=True)
            elif len(alt_nt) < len(wt_nt):  # deletion
                logging.info('deletion found at pos %d', pos)
                # indel = zip_longest(wt_nt, alt_nt, fillvalue='-')
                # for i, b in enumerate(indel):
                #     pm_here = {'wt': b[0], 'pos': pos + i, 'mut': b[1], 'freq': float(alt_f)}
                #     vcf_mutations = vcf_mutations.append(pm_here, ignore_index=True)
            pm_here = {'wt': wt_nt, 'pos': pos, 'mut': alt_nt, 'freq': float(alt_f)}
            vcf_mutations.append(pm_here)
    vcf_mutations = pd.DataFrame(vcf_mutations)

    vcf_mutations['pos'] = vcf_mutations['pos'].astype(int)
    return vcf_mutations


def merge_mutations(cns_muts, vcf_muts):
    """Merge reference bases with mutations found in a vcf file.

    Bases found on reference are passed as a DataFrame with pos/mut columns
    and frequency is intended 100%. Mutations found in vcf have a frequency
    column, so that when a mutation is found, its frequency is subtracted
    from the reference base.
    e.g. at position pos=123 consensus has T and vcf reports a T->C at 25%.
    The final table will have
    pos    mut    freq
    123     T     0.75
    123     C     0.25
    """
    mutated_pos = set(vcf_muts.pos)
    unmutated_pos = set(cns_muts.pos) - mutated_pos
    # save the unmutated positions
    unmutated = cns_muts[cns_muts.pos.isin(unmutated_pos)]
    # now add the mutated positions
    save_pos = []
    save_freq = []
    save_mut = []
    for pos in sorted(mutated_pos):
        variants_here = vcf_muts[vcf_muts.pos == pos]
        cumulative = variants_here.freq.sum()
        # save first what is found on consensus minus frequency of all vcf mutations
        save_pos.append(pos)
        save_freq.append(1. - cumulative)
        save_mut.append(cns_muts[cns_muts.pos == pos].mut.tolist()[0])
        # now add all mutations
        save_pos.extend(variants_here.pos.tolist())
        save_freq.extend(variants_here.freq.tolist())
        save_mut.extend(variants_here.mut.tolist())

    mutated = pd.DataFrame({'pos': save_pos, 'freq': save_freq, 'mut': save_mut})
    merged = pd.concat([unmutated, mutated])
    merged = merged[merged.freq > 0]
    merged = merged.sort_values(by=['pos', 'freq'], ascending=[True, False])
    return merged


def phase_mutations(bam_file, start, end):
    """Phase mutations affecting the same codon.

    If mutations appear on the same codon, it checks on the reads whether they
    really occur together.
    """
    from collections import Counter
    logging.debug('looking for mutations on the same codon to be phased')
    haps = []
    # coverage = 0
    cml = shlex.split(
        'samtools view %s sample_consensus:%d-%d' % (bam_file, start, end))
    proc = subprocess.Popen(cml, stdout=subprocess.PIPE,
                            universal_newlines=True)
    with proc.stdout as handle:
        for l in handle:
            lsp = l.split('\t')
            pos, cigar, read = int(lsp[3]), lsp[5], lsp[9]
            if 'D' in cigar or 'I' in cigar:
                continue
            # reg = read[start - pos:start - pos + 3]
            haps.append(read[start - pos:start - pos + 3])
            # coverage += 1
    logging.debug('position:%d', start)
    cnt = Counter(haps)
    coverage = sum(cnt.values())
    freq_here = [float(count) / coverage for item, count in cnt.items()]
    items_here = [item for item, count in cnt.items()]
    haps = pd.DataFrame({'mut': items_here, 'freq': freq_here})
    haps = haps[haps.freq > HAPLO_FREQ_THRESHOLD]
    renorm = haps.freq.sum()
    haps.freq = haps.freq / renorm
    # assert abs(haps.freq.sum() - 1.0) < 1.E-4, haps
    return haps


def compute_org_mutations(aa_sequence, org_found):
    """Compute mutation profile and mapping of a protein sequence.

    This function takes a protein sequence aa_sequence, aligns it to
    reference depending on org_found (H77 for org_found=HCV,
    consensus_B for org_found=HIV) and returns a DataFrame where, on
    each row, are

    - wt: amminoacid on organism sequence,
    - pos: position on the organism sequence,
    - mut: amminoacid on aa_sequence,
    - in_pos: position on aa_sequence
    """
    import re
    if org_found == 'HIV':
        org_ref = B_pol_aa_seq
    elif org_found == 'HCV':
        org_ref = h77_aa_seq

    needle_align('asis:%s' % org_ref, 'asis:%s' % aa_sequence, 'sample_vs_wt.fasta', go=40., ge=4.)
    alhr = alignfile2dict(['sample_vs_wt.fasta'])
    alih = alhr['asis']['asis']
    alih.summary()
    if 3 * alih.mismatches > alih.ident:
        warnings.warn('Too many mismatches in parsing mutations')
        logging.warning('In parsing mutations mismatches: %d ident:%d', alih.mismatches, alih.ident)

    org_seq, in_seq = alih.seq_a.upper(), alih.seq_b.upper()
    end_org, end_in = len(org_seq.rstrip('-')), len(in_seq.rstrip('-'))
    end_pos = min(end_org, end_in)

    # these will count positions on wt, mutated sequence and matching positions
    org_pos = 0
    in_pos = 0
    i = 0
    nmuts = 0
    wt_save = []
    pos_save = []
    in_save = []
    mut_save = []
    # first pass does not parse indels
    al_pairs = list(zip(org_seq[:end_pos], in_seq[:end_pos]))
    for j, z in enumerate(al_pairs[:end_pos - 1]):
        org_pos += z[0] != '-'
        in_pos += z[1] != '-'
        i += '-' not in z
        if i == 0:
            continue
        # this is not a gap, nor the start of an indel
        if '-' not in z and '-' not in al_pairs[j + 1]:
            nmuts += 1
            wt_save.append(z[0])
            pos_save.append(org_pos)
            in_save.append(in_pos)
            mut_save.append(z[1])

    try:  # treat last position explicitely
        if org_seq[end_pos] != in_seq[end_pos]:
            nmuts += 1
            wt_save.append(org_seq[end_pos])
            pos_save.append(org_pos + 1)
            in_save.append(in_pos + 1)
            mut_save.append(in_seq[end_pos])
    except IndexError:
        # if sequences end together with an ammino acid end_pos is beyond the last position
        logging.debug('sequences do not end together')

    logging.info('%d mutated sites found on sample consensus', nmuts)
    # parse insertions (i.e. gaps in org_seq)
    match = re.finditer(r"\w-+", org_seq[:end_pos])
    if match:
        match_list = list(match)
        logging.info('%d insertions found on sample consensus', len(match_list))
        for m in match_list:
            logging.debug('start:%d end:%d insertion: %s',
                          m.start(), m.end(), m.group(0))
            # pos of base before gap start, 1-based, computed on organism reference
            ins_start = len(org_seq[:m.start() + 1].replace('-', ''))
            # pos of base before gap starts, 1-based, computed on input sequence
            ins_start_in = len(in_seq[:m.start() + 1].replace('-', ''))
            wt = org_seq[m.start()]
            assert wt == org_ref[ins_start - 1]
            mut = in_seq[m.start():m.end()]
            assert '-' not in mut
            wt_save.append(wt)
            pos_save.append(ins_start)
            in_save.append(ins_start_in)
            mut_save.append(mut)

    # parse deletions (i.e. gaps in in_seq)
    match = re.finditer(r"\w-+", in_seq[:end_pos])
    if match:
        match_list = list(match)
        logging.info('%d deletions found on sample consensus', len(match_list))
        for m in match_list:
            logging.debug('start:%d end:%d deletion: %s',
                          m.start(), m.end(), m.group(0))
            # pos of base before gap_start, 1-based, computed on organism reference
            del_start = len(org_seq[:m.start() + 1].replace('-', ''))
            # pos of base before gap_start, 1-based, computed on input sequence
            del_start_in = len(in_seq[:m.start() + 1].replace('-', ''))
            wt = org_seq[m.start():m.end()]
            assert wt[0] == org_ref[del_start - 1]
            mut = m.group(0).rstrip('-')
            wt_save.append(wt)
            pos_save.append(del_start)
            in_save.append(del_start_in)
            mut_save.append(mut)

    here_muts = pd.DataFrame({'pos': pos_save, 'wt': wt_save, 'mut': mut_save, 'in_pos': in_save})
    here_muts['pos'] = here_muts['pos'].astype(int)
    return here_muts


def sequence_2_df(seq_in):
    """Return a DataFrame with pos (1-based), base and frequency (1.0)."""
    pos = [i + 1 for i in range(len(seq_in))]
    df_out = pd.DataFrame({'pos': pos, 'mut': list(seq_in), 'freq': 1.0})
    df_out['pos'] = df_out['pos'].astype(int)
    df_out['freq'] = df_out['freq'].astype(float)
    return df_out


def df_2_sequence(df_in):
    """Take a DataFrame with positions and nucleotides and returns a sequence.

    If a frequency column exists, select the positions with maximum frequency.
    """
    # select row with max frequency for each position
    if 'freq' in df_in.columns:
        max_freq_idx = df_in.groupby('pos')['freq'].transform("max") == df_in['freq']
        df_in = df_in[max_freq_idx]
    # this is needed for positions where frequency is 0.5/0.5
    df_in = df_in.groupby('pos').first().reset_index()
    df_in = df_in.sort_values(by='pos', ascending=True)
    # extract first nt in case of insertions
    df_in['single_mut'] = df_in.apply(lambda row: row['mut'][0], axis=1)
    return ''.join(df_in.single_mut.tolist())


def df_2_ambiguous_sequence(df_in, cov_df=None):
    """Take a DataFrame with positions, nucleotides and frequencies and returns a sequence.

    If the frequency is above a certain threshold, write the sequence with wobble bases.
    """
    assert 'freq' in df_in.columns
    # select calls with freq > 15%
    df_in = df_in[df_in['freq'] >= AMBIGUITY_threshold]
    # aggregate calls for the same position
    all_nt = df_in.groupby(['pos']).agg({'mut': lambda x: ''.join(sorted(x))})
    # create a columng of ambiguous bases * so insertion is not considered because only the first nc is returned
    all_nt['ambi'] = all_nt.apply(lambda row: d2a.get(row['mut'], row['mut'][0]), axis=1)
    all_nt.reset_index(inplace=True)
    if not cov_df is None:
        full_df = pd.merge(all_nt, cov_df, on='pos', how='left')
        full_df.loc[full_df['coverage'] < coverage_threshold, 'ambi'] = 'N'
    return ''.join(full_df.ambi.tolist())


def get_coverage(sorted_bam_file):
    """Run samtools depth and read the result in a pandas DataFrame."""
    cml = shlex.split('samtools depth -a %s' % sorted_bam_file)
    proc = subprocess.Popen(cml, stdout=subprocess.PIPE)
    cov_df = pd.read_csv(proc.stdout, sep='\t', header=None, names=['ref', 'pos', 'coverage'])
    return cov_df


def add_coverage_2_merged(merged, coverage_df):
    """Add the coverage information to the merged dataframe that includes 
    reference bases and mutations.
    
    This returns merged_w_cov dataframe with freq, mu, pos, coverage columns.
    """
    merged_w_cov = pd.merge(merged, coverage_df[['pos', 'coverage']], on='pos', how='left')    
    return merged_w_cov

def nt_freq_2_aa_freq(df_nt, frame, bam_file=None):
    """Compute aminoacid mutations from a DataFrame of nucleotide mutations.

    Nucleotides passed are expected to be on the same codon.
    Example::

        pos  mut    freq    --->    pos   codon  freq    --->    pos   aa   freq
        1    T      1.0              1     TTA    0.9             1     L    0.9
        2    T      0.9              1     TCA    0.1             1     S    0.9
        2    C      0.1
        3    A      1.0

    ``bam_file`` is used to phase mutations found on two positions.
    """
    from itertools import product
    df_nt = df_nt.sort_values(by=['pos', 'freq'], ascending=[True, False])
    min_pos, max_pos = min(df_nt.pos), max(df_nt.pos)
    assert max_pos - min_pos == 2, df_nt
    codon_number = 1 + int((min_pos - frame) / 3)
    assert (min_pos - frame) % 3 == 0
    low_freq_positions = set(df_nt[df_nt.freq < 1.0].pos.tolist())

    if low_freq_positions == set([]):
        logging.debug('no mutated position in codon %d', codon_number)
        try:
            cod = ''.join(df_nt.mut.tolist())
            aa = translation_table[cod]
            return [codon_number], [aa], [1.0], [cod]
        except KeyError:
            logging.warning('frameshift mutations around pos:%d', min_pos)
            return [], [], [], []
    elif len(low_freq_positions) == 1:
        logging.debug('one mutated position in codon %d', codon_number)
        freqs = df_nt[df_nt.freq < 1.0].freq.tolist()
        # build all combinations of nucleotides
        a = [df_nt[df_nt.pos == min_pos + i].mut.tolist() for i in range(3)]
        combinations = [''.join(comb) for comb in product(a[0], a[1], a[2])]

        aas = []
        save_combs = []
        for i, x in enumerate(combinations):
            if '-' in x or len(x) % 3:
                logging.warning('frameshift mutation at freq %f', freqs[i])
                aas.extend('*')
                save_combs.extend(x)
                continue
            if len(x) == 3:
                aas.extend(translation_table[x])
                save_combs.extend(x)
            else:
                nt_split = [x[i: i + 3] for i in range(0, len(x), 3)]
                aa_split = [translation_table[nth] for nth in nt_split]
                aas.extend([''.join(aa_split)])
                save_combs.extend([''.join(nt_split)])
        assert len(aas) == len(freqs), '%s - %s - %s' % (combinations, aas, freqs)
        return [codon_number] * len(freqs), aas, freqs, combinations
    else:
        logging.info('phasing needed in codon %d', codon_number)
        haps = phase_mutations(bam_file, min_pos, max_pos)
        aas = [translation_table.get(h, '*') for h in haps.mut.tolist()]
        return [codon_number] * haps.shape[0], aas, haps.freq.tolist(), haps.mut.tolist()


def gene_name_pos(df_in):
    """Return the gene according to position on the H77 polyprotein/consensus_B sequence."""
    if df_in.organism == 'HCV':
        return h77_map[df_in.pos]
    elif df_in.organism == 'HIV':
        return consensus_B_map[df_in.pos]


def main(vcf_file='hq_2_cns_final_recal.vcf', ref_file='cns_final.fasta', bam_file='hq_2_cns_final_recal.bam',
         organism='HCV'):
    """What the main does."""
    if vcf_file is None:
        for fname in ['merged_muts_drm_annotated.csv', 'cns_max_freq.fasta']:
            fhandle = open(fname, 'a')
            fhandle.close()
        return
    # ref_file is the sample consensus
    ref_nt = list(SeqIO.parse(ref_file, 'fasta'))[0]
    frame_ref, aa_framed_ref = find_frame(ref_nt.seq)
    del aa_framed_ref
    logging.info('sample consensus reference file:%s\tframe:%d', ref_file, frame_ref)
    logging.info('parsing bases on consensus reference, 1-based')
    cns_variants_nt = sequence_2_df(ref_nt)
    # cns_variants_nt.to_csv('cns_variants_nt.csv', index=False, float_format='%2.1f')

    logging.info('parsing bases in vcf file')
    vcf_stem = os.path.splitext(vcf_file)[0]
    vcf_mutations = parsevar('%s.vcf' % vcf_stem, ref_nt.seq)
    # vcf_mutations.to_csv('vcf_mutations_nt.csv', index=False, float_format='%6.4f')
    logging.info('apply mutations in vcf to sample reference')
    merged = merge_mutations(cns_variants_nt, vcf_mutations)
    coverage_df = get_coverage(bam_file)
    merged_with_coverage = add_coverage_2_merged(merged,coverage_df)
    # merged_with_coverage contains four columns: variant, frequency, position of the sample consensus and coverage
    merged_with_coverage.to_csv('merged_mutations_nt.csv', index=False, float_format='%6.4f')  
    
    logging.info('save consensus sequence with wobbles to file')
    ambi_cns = Seq(df_2_ambiguous_sequence(merged, coverage_df))
    assert len(ambi_cns) == len(ref_nt), '%d - %d' % (len(ambi_cns), len(ref_nt))
    SeqIO.write(SeqRecord(ambi_cns, id='ambiguous_consensus', description=''), 'cns_ambiguous.fasta', 'fasta')

    logging.info('save max frequency sequence to file')
    # generally slightly different from consensus reference
    max_freq_cns = Seq(df_2_sequence(merged))
    assert len(max_freq_cns) == len(ref_nt), '%d - %d' % (len(max_freq_cns), len(ref_nt))
    SeqIO.write(SeqRecord(max_freq_cns, id='max_freq_cons', description=''), 'cns_max_freq.fasta', 'fasta')
    logging.info('translate max frequency sequence')
    frame_max, aa_framed_max = find_frame(max_freq_cns)
    assert frame_max == frame_ref

    # cns_variants_aa = sequence_2_df(aa_framed)
    # cns_variants_aa.to_csv('max_freq_variants_aa.csv', index=False, float_format='%2.1f')
    logging.info('align max freq to H77/consensus_B and save mutation list')
    # max_freq_muts_aa contains two position columns:
    # - in_pos is the position on the sample consensus,
    # - pos is the position on H77/consensus_B.
    # This will be useful also later to report minority variants on H77/consensus_B
    max_freq_muts_aa = compute_org_mutations(aa_framed_max, org_found=organism)
    max_freq_muts_aa.to_csv('max_freq_muts_aa.csv', index=False, float_format='%6.4f')

    logging.info('translate minority variants going over one codon at a time')
    ref_len = len(max_freq_cns)
    last_codon_pos = ref_len - (ref_len - frame_max) % 3  # prevents the last translation from failing
    codon_positions = [(i, i + 1, i + 2) for i in range(frame_max, last_codon_pos, 3)]
    a_pos_save = []
    aas_save = []
    freqs_save = []
    nt_save = []
    nt_info_save = []
    #nt_freq_save = []
    for c in codon_positions:
        codon_muts = merged[merged.pos.isin(c)]  # extracts just the affected positions
        a_pos, aas, freqs, nts = nt_freq_2_aa_freq(codon_muts, frame_max, bam_file)
        a_pos_save.extend(a_pos)
        aas_save.extend(aas)
        freqs_save.extend(freqs)
        nt_save.extend(nts)
        #nt_freq_save.append(nt_freq)
        for i in a_pos:
            codon_muts_copy = codon_muts.copy(deep=True)
            nt_info_save.append(codon_muts_copy)
    #all_muts_aa = pd.DataFrame({'freq': freqs_save, 'in_pos': a_pos_save, 'mut': aas_save, 'nts': nt_save})
    all_muts_aa = pd.DataFrame({'freq': freqs_save, 'in_pos': a_pos_save, 'mut': aas_save, 'nts': nt_save, 'nt_info':nt_info_save })
    # join with max_freq_muts_aa on in_pos (sample consensus positions) to obtain positions on H77/consensus_B
    all_muts_aa_full = pd.merge(max_freq_muts_aa, all_muts_aa, on='in_pos')
    
    all_muts_nt_duplicate = all_muts_aa_full.loc[:,['pos','nt_info']].copy(deep=True)
    all_muts_nt_nested = all_muts_nt_duplicate.drop_duplicates(subset= 'pos')
    ref_start_pos = all_muts_nt_nested['pos'].iloc[0]
    all_muts_nt_list = all_muts_nt_nested['nt_info'].tolist()
    all_muts_nt = pd.concat(all_muts_nt_list, ignore_index=True)
    #all_muts_nt['pos'] = all_muts_nt['pos'] - all_muts_nt['pos'][0] + 1
    delta = (all_muts_nt['pos'][0] - (3 * (ref_start_pos - 1) + 1))
    all_muts_nt['pos'] = all_muts_nt['pos'] - delta
    all_muts_nt.to_csv('mutations_nt_pos_ref_aa.csv', index=False, float_format='%6.4f')
    
    assert (all_muts_nt['pos'] > 0).all() 
    if organism == 'HIV':
        #assert all_muts_nt['pos'].iloc[-1] == max_freq_muts_aa['pos'].iloc[-1] * 3
        if all_muts_nt['pos'].iloc[-1] != max_freq_muts_aa['pos'].iloc[-1] * 3:
            logging.debug("The last position of mutations_nt_pos_ref_aa.csv file did not match with the last position of max_freq_muts_aa.csv file. ")
    
    all_muts_aa_full.drop(['mut_x', 'in_pos','nt_info'], axis=1, inplace=True)
    # all_muts_aa_full.drop(['mut_x'], axis=1, inplace=True)
    all_muts_aa_full.rename(columns={'mut_y': 'mut'}, inplace=True)
    all_muts_aa_full.to_csv('intermediate.csv', index=False, float_format='%6.4f')
    all_muts_aa_codon_full_csv = all_muts_aa_full.copy(deep=True)
    all_muts_aa_full = all_muts_aa_full.drop(['nts'], axis=1)
    # sum over synonymous mutations. Since deletions look like synonymous mutations, 
    # they will be removed!
    all_muts_aa_full = all_muts_aa_full.groupby(['pos', 'wt', 'mut']).sum()
    all_muts_aa_codon_full_csv = all_muts_aa_codon_full_csv.groupby(['pos', 'wt', 'mut', 'nts']).sum()
    all_muts_aa_full.reset_index(inplace=True)
    all_muts_aa_codon_full_csv.reset_index(inplace=True)
    all_muts_aa_full['organism'] = organism
    all_muts_aa_codon_full_csv['organism'] = organism
    all_muts_aa_full['gene'], all_muts_aa_full['gene_pos'] = zip(*all_muts_aa_full.apply(gene_name_pos, axis=1))
    all_muts_aa_codon_full_csv['gene'], all_muts_aa_codon_full_csv['gene_pos'] = zip(*all_muts_aa_codon_full_csv.apply(gene_name_pos, axis=1))
    all_muts_aa_codon_full_csv.to_csv('all_muts_aa_codon_full.csv', index=False, float_format='%6.4f')
    # keep only the real mutations
    real_muts = all_muts_aa_full[all_muts_aa_full.wt != all_muts_aa_full.mut]
    real_muts = real_muts.sort_values(by=['pos', 'freq'], ascending=[True, False])
    real_muts = real_muts.drop(['organism', 'pos'], axis=1)
    real_muts = real_muts[['gene', 'wt', 'gene_pos', 'mut', 'freq']]
    real_muts.rename(columns={'gene_pos': 'pos'}, inplace=True)
    real_muts.to_csv('final.csv', index=False, float_format='%6.4f')


if __name__ == '__main__':
    main('hq_2_cns_final_recal.vcf', 'cns_final.fasta',
         'hq_2_cns_final_recal.bam', sys.argv[1])
