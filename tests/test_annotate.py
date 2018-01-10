#!/usr/bin/env python3
"""Simple test facility"""
import os

import pytest
from Bio.Seq import Seq
from Bio import SeqIO

from src.minvar.annotate import (find_frame, merge_mutations, parsevar, compute_org_mutations, nt_freq_2_aa_freq,
    df_2_sequence)

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


seq_h77_short = Seq('TTCTTTCACGACCCTGCCAGCCTTGTCCACCGGCCTCATCCACCTCCACCAGAACATTGTGGACGTGCAG')
# two mismatches
seq_h77_mism = Seq('TTCTTTCACGACCCTGCCCGCCTTGTCCACCGACCTCATCCACCTCCACCAGAACATTGTGGACGTGCAG')
# one 10 bp long insertion
seq_h77_ins = Seq('TTCTTTCACGACCCTAAAAAAAAAAGCCAGCCTTGTCCACCGGCCTCATCCACCTCCACCAGAACATTGTGGACGTGCAG')
# one short deletion, one mismatch
seq_h77_del = Seq('TTCCACGATGCCAGCCTTGTCCACCGGCCTCATCCACCTCCACCAGAACATTGTGGACGTGCAG')

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
    frame, aa_framed = find_frame(ref_1)
    assert frame == 3
    # assert len(nt_framed) % 3 == 0
    # assert str(nt_framed)[-1] == 'T'
    assert aa_framed == 'IYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFL'

@pytest.mark.skip(reason='Rewriting this')
def test_merge_mutations():
    import pandas as pd
    '''
    A, 1011, G
    A, 1050, T
    '''
    vm = pd.DataFrame([['T', 1, 'A', 0.3],    # include new mutation
                       ['G', 1011, 'T', 0.2],   # already mutated in sample cons
                       ['T', 1050, 'A', 0.1]],  # mutates back to cons B
                      columns=['wt', 'pos', 'mut', 'freq'])
    m_out = merge_mutations(vm, vm)
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


def test_org_mutations():
    short_aa = Seq('GACECPGRSRRPCTMSTNPKPQ')
    muts = compute_org_mutations(short_aa, 'HCV')
    os.remove('sample_vs_wt.fasta')
    real_muts = muts[muts.mut != muts.wt]
    assert real_muts.shape[0] <= 1

    short_aa = \
     Seq('PGPFLDKPAQCPEIWACPRKTASRVVLGRERPCGTA**GACECPGRSRRPCTMSTNPKPQRKTKRNTNRRPMDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARQPEGRSWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPNDPRRRSRNLGKVIDTLTC')
    muts = compute_org_mutations(short_aa, 'HCV')
    os.remove('sample_vs_wt.fasta')
    real_muts = muts[muts.mut != muts.wt]
    assert real_muts.shape[0] <= 5

    long_aa_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long_aa.fasta')
    long_aa = list(SeqIO.parse(long_aa_file, 'fasta'))[0].seq
    muts = compute_org_mutations(long_aa, 'HCV')
    os.remove('sample_vs_wt.fasta')
    real_muts = muts[muts.mut != muts.wt]
    assert real_muts.shape[0] <= 650


def test_nt_freq_2_aa_freq():
    import pandas as pd
    # all frames, no variants
    before = 13  # codons before this
    for frame in [1, 2, 3]:
        df = pd.DataFrame({
            'pos': [frame + 3 * before, 1 + frame + 3 * before, 2 + frame + 3 * before],
            'mut': ['A', 'A', 'A'],
            'freq': [1.0, 1.0, 1.0]
        })
        codon_pos, aas, freqs, haps = nt_freq_2_aa_freq(df, frame)
        assert codon_pos[0] == before + 1
        assert aas[0] == 'K'
        assert freqs[0] == 1.0
        assert haps[0] == 'AAA'


def test_df_2_sequence():
    import pandas as pd
    # standard case
    m = pd.DataFrame({
        'pos': [1, 2, 2, 3],
        'freq': [1.0, 0.95, 0.05, 1.0],
        'mut': ['A', 'C', 'A', 'T']
        })
    assert df_2_sequence(m) == 'ACT'
    # what happens when frequency is equally split?
    m2 = pd.DataFrame({
        'pos': [1, 2, 2, 3],
        'freq': [1.0, 0.5, 0.5, 1.0],
        'mut': ['A', 'C', 'A', 'T']
        })
    assert df_2_sequence(m2) == 'ACT'
