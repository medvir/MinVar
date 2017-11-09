#!/usr/bin/env python3
"""Interprets the vcf and translates it into amminoacid mutations."""

import logging
import os
import shlex
import subprocess
import sys
import warnings

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

dn_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if __name__ == '__main__':
    if __package__ is None:
        os.sys.path.insert(1, dn_dir)
        mod = __import__('minvar')
        sys.modules["minvar"] = mod
        from common import h77_map, consensus_B_map
        from Alignment import needle_align, alignfile2dict
else:
    from .common import h77_map, consensus_B_map
    from .Alignment import needle_align, alignfile2dict


RAW_DEPTH_THRESHOLD = 50
MIN_FRACTION = 0.015
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
    to_trim = (len(refseq) - best_frame + 1) % 3

    if to_trim:
        return best_frame, Seq(refseq[best_frame - 1:-to_trim]),\
            aa_seq[best_frame].strip('*')

    return best_frame, Seq(refseq[best_frame - 1:]),\
        aa_seq[best_frame].strip('*')


def parse_cons_mutations(input_seq, org_found):
    """Return mutations of consensus w.r.t. reference.

    It uses consensus_B for HIV, subtypre_ref.fasta for HCV. 1-based
    coordinate on input sequence because we need to phase minority calls
    later.
    """
    import re

    logging.info('parsing consensus mutations for organism %s', org_found)
    if org_found == 'HIV':
        org_ref = str(B_pol_nt_seq)
    elif org_found == 'HCV':
        org_ref = str(h77_nt_seq)
    # adapt penalties for testing with short sequences
    go = 40.0 if len(input_seq) > 100 else 10.0
    ge = 5.0 if len(input_seq) > 100 else 0.5
    needle_align('asis:%s' % org_ref,
                 'asis:%s' % str(input_seq), 'parse.tmp',
                 go=go, ge=ge)
    alhr = alignfile2dict(['parse.tmp'])
    alih = alhr['asis']['asis']
    alih.summary()
    if 3 * alih.mismatches > alih.ident:
        warnings.warn('Too many mismatches in parsing mutations')
    else:
        os.remove('parse.tmp')

    # the naming of variables here below (B_seq, end_B, ...) stems from the
    # fact that, originally, only HIV consensus_B was used
    org_seq, in_seq = alih.seq_a.upper(), alih.seq_b.upper()

    end_org, end_in = len(org_seq.rstrip('-')), len(in_seq.rstrip('-'))
    end_pos = min(end_org, end_in)

    here_muts = pd.DataFrame(columns=['wt', 'pos', 'mut'])
    # these will count positions on wt, mutated sequence and matching positions
    org_pos = 0
    in_pos = 0
    i = 0
    nmuts = 0
    freq_here = 1.0
    # first pass does not parse indels
    for z in list(zip(org_seq, in_seq))[:end_pos]:
        org_pos += z[0] != '-'
        in_pos += z[1] != '-'
        i += '-' not in z
        if i == 0:
            continue
        if z[0] != z[1] and '-' not in z:
            nmuts += 1
            dict_here = {'wt': z[0], 'pos': org_pos, 'mut': z[1],
                         'freq': freq_here}
            here_muts = here_muts.append(dict_here, ignore_index=True)

    # parse insertions (i.e. gaps in org_seq)
    match = re.finditer(r"\w-+", org_seq[:end_pos])
    if match:
        match_list = list(match)
        logging.info('%d insertions found on sample consensus', len(match_list))
        for m in match_list:
            logging.debug('start:%d end:%d insertion: %s',
                          m.start(), m.end(), m.group(0))
            # pos of base before gap start, 1-based, computed on reference
            ins_start = len(org_seq[:m.start() + 1].replace('-', ''))
            wt = org_seq[m.start()]
            assert wt == org_ref[ins_start - 1]
            mut = in_seq[m.start():m.end()]
            assert '-' not in mut
            dict_here = {'wt': wt, 'pos': ins_start, 'mut': mut,
                         'freq': freq_here}
            here_muts = here_muts.append(dict_here, ignore_index=True)

    # parse deletions (i.e. gaps in in_seq)
    match = re.finditer(r"\w-+", in_seq[:end_pos])
    if match:
        match_list = list(match)
        logging.info('%d deletions found on sample consensus', len(match_list))
        for m in match_list:
            logging.debug('start:%d end:%d deletion: %s',
                          m.start(), m.end(), m.group(0))
            # pos of base before gap_start, 1-based
            del_start = len(org_seq[:m.start() + 1].replace('-', ''))
            wt = org_seq[m.start():m.end()]
            assert wt[0] == org_ref[del_start - 1]
            mut = m.group(0).rstrip('-')
            dict_here = {'wt': wt, 'pos': del_start, 'mut': mut,
                         'freq': freq_here}
            here_muts = here_muts.append(dict_here, ignore_index=True)

    logging.info('%d mutated sites found on sample consensus', nmuts)
    here_muts['pos'] = here_muts['pos'].astype(int)

    return here_muts


def parsevar(vcf_file, ref_seq):
    """Parse mutations in vcf file.

    This also checks the effect on the translation and
    returns them together with frequencies.
    """
    import re
    from itertools import zip_longest
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

    vcf_mutations = pd.DataFrame(columns=['wt', 'pos', 'mut', 'freq'])
    # vcf_mutations = vcf_mutations.set_index(['gene', 'pos', 'mut'])
    # now parses the variants
    for vc_line in vc_lines:
        lsp = vc_line.strip().split('\t')
        if float(lsp[5]) < 100:
            logging.warning('Quality filter on position %s: qual=%s',
                            lsp[1], lsp[5])
            continue

        # assert lsp[0] == ref_seq.id, 'sequence ID not matching'
        pos = int(lsp[1])
        wt_nt, alt = lsp[3:5]
        assert ref_nt[pos - 1] == wt_nt[0], 'ref shift: %s not %s' % \
            (ref_nt[pos - 1], wt_nt)

        infos = {}
        for a in lsp[7].split(';'):
            asp = a.split('=')
            try:
                infos[asp[0]] = asp[1]
            except IndexError:
                infos[asp[0]] = None

        for alt_nt, alt_f in zip(alt.split(','), infos['AF'].split(',')):
            if len(alt_nt) != len(wt_nt):  # indel
                logging.info('indel found at pos %d', pos)
                indel = zip_longest(wt_nt, alt_nt, fillvalue='-')
                for i, b in enumerate(indel):
                    if b[0] != b[1]:
                        pm_here = {'wt': b[0], 'pos': pos + i, 'mut': b[1],
                                   'freq': float(alt_f)}
                        vcf_mutations = vcf_mutations.append(
                            pm_here, ignore_index=True)
            else:
                pm_here = {
                    'wt': wt_nt, 'pos': pos, 'mut': alt_nt, 'freq': float(alt_f)
                    }
                vcf_mutations = vcf_mutations.append(pm_here, ignore_index=True)
    vcf_mutations['pos'] = vcf_mutations['pos'].astype(int)
    return vcf_mutations


def merge_mutations(cons_muts, vcf_muts):
    """Merge mutations on samples consensus with those in vcf file.

    Mutations found on sample consensus w.r.t. organism consensus and
    minority ones found w.r.t. sample consensus must be merged to represent
    true state of variants w.r.t. consensus B. For all, position is on sample
    consensus

    e.g. let's consider
    positions    1 2 3 4 5 6 7
    org. cons    A C T G A T T
                 | |   | | |
    sample cons  A C C G T T A
                 :   :   :
    vcf muts     C   G   A

    and let's assume that both vcf variants have frequency 20%.

    In position 1 we have sample consensus conserved and a variant in vcf.
    We only need to report
    A, 1, C, 0.2

    In position 3 we have both in sample consensus and vcf a variation.
    We need to report both
    T, 3, C, 0.8 (i.e. 100% - 20%)
    and
    T, 3, G, 0.2

    In position 5 the mutation in vcf is the consensus_B wild type, so we only
    need to correct the frequency of cns and report
    A, 5, T, 0.8

    In position 7 we only have a mutation on sample consensus, so we only
    report
    T, 7, A, 1.0
    """
    case = {}
    # pos can be index on cns_mutations because they are unique
    cons_muts.set_index(['pos'], inplace=True, verify_integrity=True)
    cc = pd.DataFrame(columns=['wt', 'pos', 'mut', 'freq'])

    for k, v in vcf_muts.iterrows():
        pos = v.pos
        cns_here = cons_muts.loc[pos]

        # sample consensus is used as a reference in variant calling
        if cns_here.mut != '-' and v.wt != '-':
            assert cns_here.mut == v.wt, '%d %s %s' % (pos, cns_here, v)
        # change the reference of vcf mutations to what is found in consensus B
        v.wt = cns_here.wt

        # like pos 1 in docstring example
        if cons_muts.loc[pos, 'freq'] == 0.0:
            cc = cc.append(v)
            case[1] = case.get(1, 0) + 1
        # like pos 3 in docstring example
        elif cons_muts.loc[pos, 'freq'] > 0.0 and v.mut != cns_here.wt:
            # correct frequency of mutation in sample consensus
            cons_muts.loc[pos, 'freq'] -= v['freq']
            cc = cc.append(v)
            case[3] = case.get(3, 0) + 1
        # like pos 5 in docstring example
        elif cons_muts.loc[pos, 'freq'] > 0.0 and v.mut == cns_here.wt:
            cons_muts.loc[pos, 'freq'] -= v['freq']
            case[5] = case.get(5, 0) + 1
    logging.info('reporting mutations merging results')
    logging.info('case: %s', str(case))
    cons_muts = cons_muts.reset_index()
    cc = pd.concat([cons_muts, cc], ignore_index=True)
    cc = cc[cc.freq > 0.0]
    cc.pos = cc.pos.astype(int)
    cc = cc.sort_values(by=['pos', 'freq'], ascending=[True, False])
    return cc


def phase_mutations(muts, frame, bam_file):
    """Phase mutations affecting the same codon.

    If mutations appear on the same codon, it checks on the reads whether they
    really occur together.
    """
    from operator import itemgetter

    logging.info('looking for mutations on the same codon to be phased')
    # find mutation positions appearing together in a codon, save in targets
    pm = pd.DataFrame(columns=['wt', 'pos', 'mut', 'freq'])
    # mut_pos = set(muts.pos)
    targets = []
    # list all codons
    for i in range(frame, 15000, 3):
        codon_pos = set([i, i + 1, i + 2])
        # mutated positions overlapping the codon
        mut_over = muts[muts.pos.isin(codon_pos)]

        if mut_over.empty:
            continue

        # positions where a deletion is detected
        del_muts = mut_over[mut_over.mut == '-']
        # if frameshift, nothing is saved
        if len(del_muts) < 3:
            # warnings.warn('frameshift del detected pos:%d-%d' % (i, i + 2))
            logging.debug(
                'frameshift deletion detected at pos:%d-%d', i, i + 2)
        # if in frame, save with mean frequency
        elif len(del_muts) == 3:
            freq_here = sum(del_muts.freq) / 3
            wt_codon = ''.join(del_muts.sort_values(by='pos').wt)
            print('deletion detected: %s%d-' % (wt_codon, i))
            d_here = {'wt': wt_codon, 'pos': i, 'mut': '---', 'freq': freq_here}
            pm = pm.append(d_here, ignore_index=True)
            # target is anyway appended to check the non indel mutations
            targets.append((i, i + 2))
            continue

        # positions where an insertion is detected
        ins_muts = mut_over[mut_over.wt == '-']
        # if frameshift, nothing is saved
        if len(ins_muts) < 3:
            # warnings.warn('frameshift ins detected pos:%d-%d' % (i, i + 2))
            logging.debug(
                'frameshift insertion detected at pos:%d-%d', i, i + 2)
        # if in frame, save with mean frequency
        elif len(ins_muts) == 3:
            print('insertion detected at pos:%s' % i)
            freq_here = sum(ins_muts.freq) / 3
            mut_codon = ''.join(ins_muts.sort_values(by='pos').mut)
            d_here = {
                'wt': '---', 'pos': i, 'mut': mut_codon, 'freq': freq_here}
            pm = pm.append(d_here, ignore_index=True)
            # target is anyway appended to check the non indel mutations
            targets.append((i, i + 2))
            continue

        # this is reached when no in frame-indel is detected
        if len(mut_over) > 1:
            targets.append((i, i + 2))

    logging.info('%d targets found', len(targets))
    # count three base haplotypes with samtools
    for t in targets:
        haps = {}
        coverage = 0
        cml = shlex.split(
            'samtools view %s padded_consensus:%d-%d' % (bam_file, t[0], t[1]))
        proc = subprocess.Popen(cml, stdout=subprocess.PIPE,
                                universal_newlines=True)
        with proc.stdout as handle:
            for i, l in enumerate(handle):
                lsp = l.split('\t')
                pos = int(lsp[3])
                cigar = lsp[5]
                if 'D' in cigar or 'I' in cigar:
                    continue
                read = lsp[9]
                reg = read[t[0] - pos:t[0] - pos + 3]
                haps[reg] = haps.get(reg, 0) + 1
                coverage += 1
        logging.info('position:%d', t[0])
        for k, v in sorted(haps.items(), key=itemgetter(1), reverse=True):
            freq_here = float(v) / coverage
            if freq_here > HAPLO_FREQ_THRESHOLD and len(k) > 1:
                for disp in range(3):
                    try:
                        muts = muts[muts.pos != t[0] + disp]
                    except ValueError:
                        pass
                d_here = {'wt': 'XYZ', 'pos': t[0], 'mut': k, 'freq': freq_here}
                pm = pm.append(d_here, ignore_index=True)
                logging.info('hap: %s freq:%f', k, freq_here)
        logging.info('<------>')
    muts = muts.reset_index()
    pm = pm.reset_index()
    phased_muts = pd.concat([muts, pm], ignore_index=True)
    phased_muts.pos = phased_muts.pos.astype(int)
    phased_muts = phased_muts.sort_values(
        by=['pos', 'freq'], ascending=[True, False])
    phased_muts = phased_muts.drop('index', axis=1)
    return phased_muts


def extract_hcv_orf(gen_ref_seq):
    """Extract the HCV ORF from the given sequence.

    HCV genome starts with a 341 nt 5'UTR region that is quite well
    conserved. After this region starts the long ORF where all
    proteins are encoded. This function cuts the input sequence from
    the beginning of the ORF to the end of the genome by aligning to
    the 5'UTR sequence.
    """
    needle_align(
        'asis:%s' % str(utr_5p_seq), 'asis:%s' % str(gen_ref_seq),
        'z.tmp', go=40.0, ge=5.0)
    alhr = alignfile2dict(['z.tmp'])
    os.remove('z.tmp')
    alih = alhr['asis']['asis']
    alih.summary()
    utr, ref = alih.seq_a.upper(), alih.seq_b.upper()
    end_utr = len(utr.rstrip('-'))
    cut_ref = ref[end_utr:]
    return Seq(cut_ref), end_utr


def get_align_map(ref_nt, B_nt):

    frame, framed_nt, framed_aa = find_frame(ref_nt)
    del frame
    del framed_nt
    # build map of sample consensus to DRM/RAS reference
    needle_align('asis:%s' % str(B_nt.translate()),
                 'asis:%s' % str(framed_aa), 'h.tmp',
                 go=40.0, ge=5.0)
    alhr = alignfile2dict(['h.tmp'])
    alih = alhr['asis']['asis']
    alih.summary()
    if 3 * alih.mismatches > alih.ident:
        warnings.warn('Too many mismatches in ref vs consensus B')
    else:
        os.remove('h.tmp')
    ref, cur = alih.seq_a.upper(), alih.seq_b.upper()
    end_ref, end_cur = len(ref.rstrip('-')), len(cur.rstrip('-'))
    end_pos = min(end_ref, end_cur)
    # this saves map[pos on sample reference] = pos on DRM/RAS nucleotide ref
    ref_pos = 0
    extend_ref = []
    for z in ref[:end_pos]:
        if z != '-':
            coords = [str(p) for p in range(ref_pos, ref_pos + 3)]
            ref_pos += 3
        else:
            coords = ['-', '-', '-']
        extend_ref.extend(coords)
    ref_pos = 0
    extend_cur = []
    for z in cur[:end_pos]:
        if z != '-':
            coords = [str(p) for p in range(ref_pos, ref_pos + 3)]
            ref_pos += 3
        else:
            coords = ['-', '-', '-']
        extend_cur.extend(coords)
    align_map = {}
    ref_pos = 0
    cur_pos = 0
    for z in zip(extend_cur, extend_ref):
        cur_pos += z[0] != '-'
        ref_pos += z[1] != '-'
        if ref_pos == 0:
            align_map[cur_pos] = '0-utr'
            continue
        if z[0] == '-':
            # deletion
            align_map[cur_pos] = '%d-del' % ref_pos
        elif z[1] == '-':
            # insertion
            align_map[cur_pos] = '%d-ins' % ref_pos
        else:
            align_map[cur_pos] = '%d' % ref_pos

    return align_map


def place_mutation(mutation, nt_ref, alignment_map):
    """Mutation nt_ref with mutation, knowing the alignment map.

    A few examples, assuming that mutation.pos is 4:

    ## Insertion

    1 2 3       4 5 6 7 8 9
    A C G - - - C G A C T T  nt_ref
    A C C A C C C G G C T G  sample consensus
    1 2 3 4 5 6 7 8 9 0 1 2

    This should result in:
    A C G A C C C G A C T T

    ## Deletion

    1 2 3 4 5 6 7 8 9 0 1 2
    A C G C G A C T T C T G  nt_ref
    A C C - - - C G G C T G  sample consensus
    1 2 3       4 5 6 7 8 9

    This should result in:
    A C G - - - C T T C T G

    ## Mismatch

    1 2 3 4 5 6 7 8 9
    A C G C G A C T G  nt_ref
    A C C C G G C T A  sample consensus
    1 2 3 4 5 6 7 8 9

    This should result in:
    A C G - - - C T T C T G



     """

    alt_nt = mutation.mut
    ref_pos = alignment_map[mutation.pos]


def annotate_mutations(mutations, ref, org_found):
    """Write gene, position, wild type, mutation, frequency in a data frame.

    Takes a DataFrame with nucleotide mutations and translates and
    annotates them according to coordinates in consensus B genes.
    Position on mutations is on ref (sample consensus), so we need to
    report these onto consensus B or HCV reference.

    EDIT:
    This function works on both HIV and HCV, but the naming of the
    variables reflect the fact that it was written for HIV only at the
    beginning.
    """
    anno_variants = pd.DataFrame(columns=['gene', 'pos', 'wt', 'mut', 'freq'])
    if mutations.empty:
        return anno_variants
    ref_nt = ref.seq

    # save the DRM/RAS reference
    if org_found == 'HIV':
        B_nt = B_pol_nt_seq
        B_aa = B_pol_aa_seq
        # orf_pos = 0
    elif org_found == 'HCV':
        # sub_ref_nt_seq = list(SeqIO.parse('subtype_ref.fasta','fasta'))[0].seq
        B_nt = h77_nt_seq
        B_aa = h77_aa_seq
        # orf_seq, orf_pos = extract_hcv_orf(sub_ref_nt_seq)
        # del orf_seq
        # del frame
    align_map = get_align_map(ref_nt, B_nt)

    # annotate each mutation
    for k, v in mutations.iterrows():
        del k
        # position on nt DRM/RAS reference
        mapped = align_map[v.pos]
        if mapped == '0-utr':
            continue
        try:
            ref_pos = int(mapped)
            pattern = 'match'
        except ValueError:
            ref_pos, pattern = int(mapped.split('-')[0]), mapped.split('-')[1]
        alt_nt = v.mut
        # mutate the DRM/RAS reference in this specific position
        if len(v.wt) == 1 and len(v.mut) == 1:  # v.wt can be XYZ after phasing
            print(ref_pos, pattern, mapped)
            print(B_nt[ref_pos - 2:ref_pos + 2])
            #assert B_nt[ref_pos - 1] == v.wt, \
            #    '%s is not %s' % (B_nt[ref_pos - 1], v)
        mut_nt = B_nt[:ref_pos - 1] + alt_nt + B_nt[ref_pos - 1 + len(alt_nt):]
        mut_nt = place_mutation(v, B_nt, align_map)
        # mut_nt is the whole sequence with just one mutated codon
        # B_nt is already in frame -> mut_nt must be in frame too
        # translate it all
        try:
            mut_aa = mut_nt.translate()
        except Bio.Data.CodonTable.TranslationError:
            if '---' in mut_nt:
                mut_aa = ''.join(
                    [translation_table[str(mut_nt[i:i + 3])]
                     for i in range(0, len(mut_nt), 3)])
            else:
                warnings.warn(
                    'CodonTable transl. err.: mut_nt:%s wt:%s pos:%d alt:%s' %
                    (mut_nt, v.wt, v.pos, alt_nt))
                continue
        if mut_aa == B_aa:
            # silent mutation
            continue
        assert len(mut_aa) == len(B_aa), \
            '%d is not %d' % (len(mut_aa), len(B_aa))
        # 0-based position of the mutated codon on prot DRM/RAS reference
        aa_pos = int((ref_pos - 1) / 3)
        assert aa_pos <= len(mut_aa), '%d %d %s' % (aa_pos, len(mut_aa), mut_aa)

        a_mut = []
        for i, a in enumerate(zip(B_aa, mut_aa)):
            if a[0] != a[1]:
                a_mut.append((i, a))
        assert len(a_mut) == 1, a_mut
        if aa_pos != a_mut[0][0]:
            warnings.warn('Is it %d or %d?' % (a_mut[0][0], aa_pos))
        B_codon, mut_codon = a_mut[0][1]
        if B_codon == mut_codon:
            warnings.warn('There must be a mutation here')
        B_pos = aa_pos + 1  # 1-based postion on prot DRM/RAS reference
        if org_found == 'HIV':
            gene_name, gene_pos = consensus_B_map[B_pos]
        elif org_found == 'HCV':
            gene_name, gene_pos = h77_map[B_pos]
        if gene_pos <= 0:
            continue
        # mut_tp = '%s%d%s' % (z[0], i, z[1])
        dict_here = {'gene': gene_name, 'wt': B_codon,
                     'pos': gene_pos, 'mut': mut_codon, 'freq': v.freq}
        anno_variants = anno_variants.append(dict_here, ignore_index=True)

    anno_variants['pos'] = anno_variants['pos'].astype(int)

    return anno_variants


def HIV_gene_map(pos):
    """Return the gene name and the position on the gene of a codon found."""
    if pos <= 56:
        gene_name = 'GagPolTF'
        gene_pos = pos
    elif 57 <= pos and pos <= 155:
        gene_name = 'protease'
        gene_pos = pos - 56
    elif 156 <= pos and pos <= 595:
        gene_name = 'RT'
        gene_pos = pos - 155
    elif 596 <= pos and pos <= 715:
        gene_name = 'RNase'
        gene_pos = pos - 595
    elif 716 <= pos and pos <= 1003:
        gene_name = 'integrase'
        gene_pos = pos - 715

    return gene_name, gene_pos


def main(vcf_file=None, ref_file=None, bam_file=None, organism=None):
    """What the main does."""
    # ref_file is the sample consensus
    # parse mutations already found in the sample consensus wrt to
    # organism reference (consensus B pol gene for HIV)
    ref_nt = list(SeqIO.parse(ref_file, 'fasta'))[0]
    frame, nt_framed, aa_framed = find_frame(ref_nt.seq)
    del nt_framed
    del aa_framed
    logging.info(
        'sample consensus reference file:%s\tframe:%d', ref_file, frame)
    logging.info('parsing mutations on consensus')
    cns_mutations = parse_cons_mutations(ref_nt, organism)
    # cns_mutations also has positions without mutations, used later
    # to merge those in vcf_mutations
    # only write real mutations to file
    cns_diff = cns_mutations.set_index('pos')
    cns_diff = cns_diff[cns_diff.freq == 1.0]
    cns_diff.to_csv('cns_mutations_nt.csv', index=True)

    vcf_stem = os.path.splitext(vcf_file)[0]

    # parse mutations in vcf file
    vcf_mutations = parsevar('%s.vcf' % vcf_stem, ref_nt)

    if vcf_mutations.shape[0]:
        # synonymous mutations can be merged
        vcf_mutations = vcf_mutations.groupby(
            ['pos', 'wt', 'mut']).sum().reset_index()
    # vcf_mutations.set_index(['pos', 'wt'], inplace=True)
    vcf_mutations.to_csv(
        'vcf_mutations_nt.csv', index=False, float_format='%6.4f')

    # now we have cns_mutations, retrieved from consensus sequence in
    # a fasta file and vcf_mutations, retrieved from vcf file created
    # with a variant calling method; while cns_mutations has mutations
    # of cns_final.fasta w.r.t. consensus B, vcf_mutations has minority
    # mutations w.r.t. cns_final.fasta. We need to subtract vcf
    # frequencies from cns ones the result is mutations w.r.t.
    # consensus B, position is on cns_final
    merged = merge_mutations(cns_mutations, vcf_mutations)

    # check frequencies again
    ch = merged.groupby(['pos', 'wt', 'mut']).sum()
    msk = ch.freq <= 1.0
    if not msk.all():
        warnings.warn('frequencies should be normalised')

    merged.to_csv('merged_mutations_nt.csv', index=False, float_format='%6.4f')

    # another step to phase variants that occur together on reads, kind of
    # making haplotypes, but only three nt long (one codon)
    phased = phase_mutations(merged, frame, bam_file)
    phased.to_csv('phased.csv', sep=',', float_format='%6.4f', index=False)

    # mutations can now be annotated and saved
    anno_muts = annotate_mutations(phased, ref_nt, organism)
    anno_muts = anno_muts.groupby(['gene', 'pos', 'wt', 'mut']).sum()
    anno_muts = anno_muts.reset_index()
    anno_muts = anno_muts.sort_values(
        by=['gene', 'pos', 'freq'], ascending=[True, True, False])
    anno_muts.to_csv(
        'annotated_mutations.csv', sep=',', float_format='%6.4f', index=False)

if __name__ == '__main__':
    main('hq_2_cns_final_recal.vcf', 'cns_final.fasta',
         'hq_2_cns_final_recal.bam', sys.argv[1])
