#!/usr/local/bin/python3.4
import pysam
import warnings

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

    return best_frame, aa_seq[best_frame]


def my_fetch_callback( alignment ):
    print(str(alignment))

def my_pileup_callback( pileups ):
    print(str(pileups))

class Counter:
    def __init__(self):
        self.counts = 0
    def __call__(self, alignment):
        self.counts += 1


def callvar(frame, framed_ref, bamfile):
    ''''''
    from collections import Counter
    bam = pysam.Samfile(bamfile, 'rb')
    for k, v in bam.header.items():
        print(k, v)
    assert len(bam.references) == 1
    ref_name = bam.references[0]

    indelreads = 0
    noindelreads = 0
    for x in bam.fetch(ref_name):
        # reads with indels have a different treatment
        if 'D' in x.cigarstring or 'I' in x.cigarstring:
            assert(len(x.blocks) > 1)
            indelreads += 1
        else:
            assert(len(x.blocks) == 1)
            noindelreads += 1
            totrans = x.query[frame - 1:]
            quals = x.qqual[frame - 1:]
#            translated = [totrans[i:i + 3].decode('utf-8') \
#                          for i in range(0, int(len(totrans) / 3) * 3, 3)]
            for i in range(0, int(len(totrans) / 3) * 3, 3):
                res = totrans[i:i + 3].decode('utf-8')
                qq = quals[i:i + 3]
                print(translation_table[res], end='\t')
                for q in qq:
                    print(q)#int.from_bytes(q, byteorder='big'), end='\t')
                print()
            break

    print(indelreads, noindelreads)

    bam.close()

def main(ref_file='matching_reference.fasta', bamfile='hq_2_cons_sorted.bam'):
    ''''''
    frame, framed_ref = find_frame(ref_file)
    callvar(frame, framed_ref, bamfile)

if __name__ == '__main__':
    main()