#!/usr/bin/env python3
"""Simple test facility"""
import os

import pytest
from Bio.Seq import Seq

from src.minvar.annotate import (annotate_mutations, extract_hcv_orf,
                                 find_frame, merge_mutations,
                                 parse_cons_mutations, parsevar)


subtype1a = '''GCCAGCCCCCTGATGGGGGCGACACTCCACCATGAATCACTCCCCTGTGAGGAACTACTG\
TCTTCACGCAGAAAGCGTCTAGCCATGGCGTTAGTATGAGTGTCGTGCAGCCTCCAGGAC\
CCCCCCTCCCGGGAGAGCCATAGTGGTCTGCGGAACCGGTGAGTACACCGGAATTGCCAG\
GACGACCGGGTCCTTTCTTGGATCAACCCGCTCAATGCCTGGAGATTTGGGCGTGCCCCC\
GCAAGACTGCTAGCCGAGTAGTGTTGGGTCGCGAAAGGCCTTGTGGTACTGCCTGATAGG\
GTGCTTGCGAGTGCCCCGGGAGGTCTCGTAGACCGTGCACCATGAGCACGAATCCTAAAC\
CTCAAAAAAAAAACAAACGTAACACCAACCGTCGCCCACAGGACGTCAAGTTCCCGGGTG\
GCGGTCAGATCGTTGGTGGAGTTTACTTGTTGCCGCGCAGGGGCCCTAGATTGGGTGTGC\
GCGCGACGAGAAAGACTTCCGAGCGGTCGCAACCTCGAGGTAGACGTCAGCCTATCCCCA\
AGGCTCGTCGGCCCGAGGGCAGGACCTGGGCTCAGCCCGGGTACCCTTGGCCCCTCTATG\
GCAATGAGGGCTGCGGGTGGGCGGGATGGCTCCTGTCTCCCCGTGGCTCTCGGCCTAGCT\
GGGGCCCCACAGACCCCCGGCGTAGGTCGCGCAATTTGGGTAAGGTCATCGATACCCTTA\
CGTGCGGCTTCGCCGACCTCATGGGGTACATACCGCTCGTCGGCGCCCCTCTTGGAGGCG\
CTGCCAGGGCCCTGGCGCATGGCGTCCGGGTTCTGGAAGACGGCGTGAACTATGCAACAG\
GGAACCTTCCTGGTTGCTCTTTCTCTATCTTCCTTCTGGCCCTGCTCTCTTGCTTGACTG\
TGCCCGCTTCGGCCTACCAAGTGCGCAACTCCACGGGGCTTTACCACGTCACCAATGATT\
GCCCTAACTCGAGTATTGTGTACGAGGCGGCCGATGCCATCCTGCACACTCCGGGGTGCG\
TCCCTTGCGTTCGTGAGGGCAACGCCTCGAGGTGTTGGGTGGCGATGACCCCTACGGTGG\
CCACCAGGGATGGCAAACTCCCCGCGACGCAGCTTCGACGTCACATCGATCTGCTTGTCG\
GGAGCGCCACCCTCTGTTCGGCCCTCTACGTGGGGGACCTATGCGGGTCTGTCTTTCTTG\
TCGGCCAACTGTTCACCTTCTCTCCCAGGCGCCACTGGACGACGCAAGGTTGCAATTGCT\
CTATCTATCCCGGCCATATAACGGGTCACCGCATGGCATGGGATATGATGATGAACTGGT\
CCCCTACGACGGCGTTGGTAATGGCTCAGCTGCTCCGGATCCCACAAGCCATCTTGGACA\
TGATCGCTGGTGCTCACTGGGGAGTCCTGGCGGGCATAGCGTATTTCTCCATGGTGGGGA\
ACTGGGCGAAGGTCCTGGTAGTGCTGCTGCTATTTGCCGGCGTCGACGCGGAAACCCACG\
TCACCGGGGGAAGTGCCGGCCACACTGTGTCTGGATTTGTTAGCCTCCTCGCACCAGGCG\
CCAAGCAGAACGTCCAGCTGATCAACACCAACGGCAGTTGGCACCTCAATAGCACGGCCC\
TGAACTGCAATGATAGCCTCAACACCGGCTGGTTGGCAGGGCTTTTCTATCACCACAAGT\
TCAACTCTTCAGGCTGTCCTGAGAGGCTAGCCAGCTGCCGACCCCTTACCGATTTTGACC\
AGGGCTGGGGCCCTATCAGTTATGCCAACGGAAGCGGCCCCGACCAGCGCCCCTACTGCT\
GGCACTACCCCCCAAAACCTTGCGGTATTGTGCCCGCGAAGAGTGTGTGTGGTCCGGTAT\
ATTGCTTCACTCCCAGCCCCGTGGTGGTGGGAACGACCGACAGGTCGGGCGCGCCCACCT\
ACAGCTGGGGTGAAAATGATACGGACGTCTTCGTCCTTAACAATACCAGGCCACCGCTGG\
GCAATTGGTTCGGTTGTACCTGGATGAACTCAACTGGATTCACCAAAGTGTGCGGAGCGC\
CTCCTTGTGTCATCGGAGGGGCGGGCAACAACACCCTGCACTGCCCCACTGATTGCTTCC\
GCAAGCATCCGGACGCCACATACTCTCGGTGCGGCTCCGGTCCCTGGATCACACCCAGGT\
GCCTGGTCGACTACCCGTATAGGCTTTGGCATTATCCTTGTACCATCAACTACACCATAT\
TTAAAATCAGGATGTACGTGGGAGGGGTCGAACACAGGCTGGAAGCTGCCTGCAACTGGA\
CGCGGGGCGAACGTTGCGATCTGGAAGACAGGGACAGGTCCGAGCTCAGCCCGTTACTGC\
TGACCACTACACAGTGGCAGGTCCTCCCGTGTTCCTTCACAACCCTACCAGCCTTGTCCA\
CCGGCCTCATCCACCTCCACCAGAACATTGTGGACGTGCAGTACTTGTACGGGGTGGGGT\
CAAGCATCGCGTCCTGGGCCATTAAGTGGGAGTACGTCGTTCTCCTGTTCCTTCTGCTTG\
'''


def test_extract_1a():
    """Identify where orf is: test with genotype 1a"""
    s, c = extract_hcv_orf(Seq(subtype1a))
    sta = str(s.translate())
    assert sta[0] == 'M'


# no mismatches
seq_h77_short = Seq('TTCTTTCACGACCCTGCCAGCCTTGTCCACCGGCCTCATCCACCTCCACCAGAACATTGTGGACGTGCAG')
# two mismatches
seq_h77_mism = Seq('TTCTTTCACGACCCTGCCCGCCTTGTCCACCGACCTCATCCACCTCCACCAGAACATTGTGGACGTGCAG')
# one 10 bp long insertion
seq_h77_ins = Seq('TTCTTTCACGACCCTAAAAAAAAAAGCCAGCCTTGTCCACCGGCCTCATCCACCTCCACCAGAACATTGTGGACGTGCAG')
# one short deletion, one mismatch
seq_h77_del = Seq('TTCCACGATGCCAGCCTTGTCCACCGGCCTCATCCACCTCCACCAGAACATTGTGGACGTGCAG')

@pytest.mark.skip(reason="rewritten")
def test_parse_cons_mutations():

    # test_no_mutations
    m = parse_cons_mutations(seq_h77_short, 'HCV')
    assert m.empty
    # test_two_mismatches
    m = parse_cons_mutations(seq_h77_mism, 'HCV')
    assert len(m) == 2
    # test_one_insertion
    m = parse_cons_mutations(seq_h77_ins, 'HCV')
    assert len(m) == 1
    # test_one_del_one_mism
    m = parse_cons_mutations(seq_h77_del, 'HCV')
    assert len(m) == 2
    assert m.iloc[0].wt == 'T'
    assert m.iloc[0].mut == 'C'


# two silent mutations w.r.t. consensus B
ref_1 = Seq('TTATCTATCAGTACATGGATGATTTGTATGTAGGATCTGACTTAGAAATTGGGCAGCATAGAACAAAAATAGAGGAACTGAGACAACATCTGTTGAGGTGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTT')
# two mutation w.r.t. consensus B, one silent
ref_2 = Seq('TTATCTATCAGTACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGACAACATCTGTTGAGGTGCGGATTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTT')
# frame is 1 here
ref_3 = Seq('TTTTTTAGGGAAGATCTGGCCTTCCCACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAAGAGACAACAACTCCCTCTCAGAAGCAGGAGCCGATAGACAAGGA')

test_vcf = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'vcf_test.vcf')
assert os.path.exists(test_vcf)


def test_parsevar():
    out = parsevar(test_vcf, ref_1)
    assert len(out) == 4, out
    out_0 = out.iloc[0]
    out_1 = out.iloc[1]

    assert out_0['wt'] == 'T', out
    assert out_0['pos'] == 4, out
    assert out_0['mut'] == 'A', out
    assert out_1['wt'] == 'G', out
    assert out_1['pos'] == 11, out
    assert out_1['mut'] == 'C', out
    assert out.iloc[2]['wt'] == 'G', out
    assert out.iloc[2]['mut'] == 'CAAA', out
    assert out.iloc[3]['wt'] == 'TAGG', out
    assert out.iloc[3]['mut'] == 'C', out


def test_find_frame():
    frame, nt_framed, aa_framed = find_frame(ref_1)
    assert frame == 3
    assert len(nt_framed) % 3 == 0
    assert str(nt_framed)[-1] == 'T'
    assert aa_framed == 'IYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFL'

@pytest.mark.skip(reason='Rewriting this')
def test_merge_mutations():
    import pandas as pd
    cm = parse_cons_mutations(ref_1, 'HIV')
    '''
    A, 1011, G
    A, 1050, T
    '''
    vm = pd.DataFrame([['T', 1, 'A', 0.3],    # include new mutation
                       ['G', 1011, 'T', 0.2],   # already mutated in sample cons
                       ['T', 1050, 'A', 0.1]],  # mutates back to cons B
                      columns=['wt', 'pos', 'mut', 'freq'])
    m_out = merge_mutations(cm, vm)
    # two mutations in consensus, three in vcf but the one in pos 101 is
    # reverting to consensus B
    assert len(m_out) == 4
    # freqs = [0.3, 0.8, 0.2, 0.9]
    # muts = ['A', 'G', 'T', 'T']
    # for i in range(4):
    #     out_here = m_out.iloc[i]
    #     with self.subTest(mut=i):
    #         self.assertEqual(out_here.freq, freqs[i])
    #         self.assertEqual(out_here.mut, muts[i])


@pytest.mark.skip(reason="check this")
def test_annotate_mutations():
    import pandas as pd
    muts = pd.DataFrame(columns=['gene', 'pos', 'wt', 'mut', 'freq'])

    # no mutation first
    out = annotate_mutations(muts, ref_3)
    assert out.empty

    # add one mutation
    muts = muts.append({'wt': 'T', 'pos': 6, 'mut': 'A', 'freq': 0.2}, ignore_index=True)
    out = annotate_mutations(muts, ref_3)
    assert len(out) == 1
    out_0 = out.iloc[0]
    assert out_0.pos == 2
    assert out_0.gene == 'GagPolTF'
    assert out_0.wt == 'F'
    assert out_0.mut == 'L'
    assert out_0.freq == 0.2
