#!/usr/local/bin/python3.4
import pysam
from Bio import SeqIO

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


def find_frame(ref):
    ''''''
    import re

    max_len = -1
    refseq = str(list(SeqIO.parse(ref, 'fasta'))[0].seq).upper()
    aa_seq = {}
    for f in [1, 2, 3]:
        nts = refseq[f-1:]
        codons = [nts[i:i+3] for i in range(0, len(nts), 3)]
        aa_seq[f] = ''.join([translation_table.get(c, '*') for c in codons])
        
        # get maximal length of ORFs
        len1 = max([len(a) for a in aa_seq[f].split('*')])
        print(len1)
        if len1 > max_len:
            best_frame = f
            
            max_len = len1
    return best_frame, aa_seq[best_frame]


def main(ref='matching_reference.fasta', bamfile='hq_2_cons_sorted.bam'):
    ''''''
    frame, framed_ref = find_frame(ref)


if __name__ == '__main__':
    main()