#!/usr/local/bin/python3.4

import sys
import subprocess
import os
import warnings
import pysam

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

RAW_DEPTH_THRESHOLD = 20
MAPPING_QUALITY_THRESHOLD = 50

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

'''
These are the fields in vcf defined by samtools mpileup.
Since using freebayes these are not needed anymore
info_descriptions = {  
	'INDEL': 'Indicates that the variant is an INDEL.',
	'IDV': 'Maximum number of reads supporting an indel',
	'IMF': 'Maximum fraction of reads supporting an indel',
    'DP': 'Raw read depth',
    'VDB': 'Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)',
	'RPB': 'Mann-Whitney U test of Read Position Bias (bigger is better)',
	'MQB': 'Mann-Whitney U test of Mapping Quality Bias (bigger is better)',
	'BQB': 'Mann-Whitney U test of Base Quality Bias (bigger is better)',
	'MQSB': 'Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)',
	'SGB': 'Segregation based metric.',
	'MQ0F': 'Fraction of MQ0 reads (smaller is better)',
	'I16': 'Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h',
	'QS': 'Auxiliary tag used for calling',
	'ICB': 'Inbreeding Coefficient Binomial test (bigger is better)',
	'HOB': 'Bias in the number of HOMs number (smaller is better)',
	'AC': 'Allele count in genotypes for each ALT allele, in the same order as listed',
	'AN': 'Total number of alleles in called genotypes',
	'DP4': 'Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases',
	'MQ': 'Average mapping quality',
}
'''
class Variant():
    '''Simple container for the called variant according to VCF specs (4.2)'''
    def __init__(self, position=0, wt=None, mut=None, chr_name=None):
        self.alt = mut  # mutatis (to)
        self.ref = wt  # mutandis (from)
        self.pos = position  # [integer]
        self.chrom = chr_name  # reference sequence [string]
        self.comment = None
        self.id = 'None'  # string
        self.qual = 'None'  # integer
        self.filter = []  # strings
        self.info = []
        self.fwd_c = 0  # number of observations on forward strand
        self.rev_c = 0  # number of observations on reverse strand
        self.coverage = 0  # total coverage 

    def update_count(self, ar_rev):
        if ar_rev:
            self.rev_c += 1
        else:
            self.fwd_c += 1

    def __str__(self):
        if self.ref:
            return '%s%d%s' % (self.ref, self.pos, self.alt)
        else:
            return '%d%s' % (self.pos, self.alt)

    def tovcf(self):
        fields = [self.chrom, str(self.pos), self.id, self.ref, self.alt,
                  str(self.qual), ';'.join(self.filter), ';'.join(self.info)]
        return '\t'.join(fields)


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


def strand_bias(freqs):
    ''''''
    return False

def filter_variant(info_field):
    '''Parses the info field in vcf file and returns false if one filter is
    triggered'''
    infos = dict(a.split('=') for a in info_field.split(';'))

    if int(infos['DP']) < RAW_DEPTH_THRESHOLD:
        return False

    freqs = [float(f) for f in infos['AF'].split(',')]
    # insert frequency of the reference at the beginning
    f0 = 1.0 - sum(freqs)
    freqs.insert(0, f0)

    return freqs

def parsevar(frame, nt_ref, aa_ref, vcf_file):
    ''''''
    import re

    # nt_ref and aa_ref are already framed
    Dn = 0
    Ds = 0
    allvars = []

    prog = re.compile('##INFO=<ID=(.*),Number.*,Description="(.*)">')

    all_muts = []
    with open(vcf_file) as h:
        for l in h:
            if l.startswith('#'):
                res = prog.match(l)
                try:
                    info_fields[res.group(1)] = res.group(2)
                    continue
                except:
                    continue
            lsp = l.split('\t')
            pos = int(lsp[1])

            assert lsp[3][0] == nt_ref[pos - frame], 'reference shift: %s' % lsp

            # pos - frame is 0-based position on the framed nt reference
            # consequently, (pos - frame) / 3 is on aa reference
            codon_pos = int((pos - frame) / 3)

            # triplet will be on positions (codon_pos * 3, (codon_pos + 1) * 3)
            wt_cod = nt_ref[codon_pos * 3: (codon_pos + 1) * 3]
            assert translation_table[wt_cod] == aa_ref[codon_pos]

            try:
                freqs = filter_variant(lsp[7])
            except:
                for k, v in info_fields.items():
                    print(k, v, file=sys.stderr)
                sys.exit('Accepted info fields are printed above')
            if not freqs:
                continue
            alt = lsp[4]

            # loop needed when there is more than one allele
            if freqs[0]:
                sr = SeqRecord(Seq(aa_ref), id='ref;freq=%f' % freqs[0],
                               description='')
                all_muts.append(sr)
            for i, alt_nt in enumerate(alt.split(',')):
                mut_nt = nt_ref[:pos - frame] + alt_nt
                mut_nt += nt_ref[pos - frame + len(alt_nt):]
                mut_seq = Seq(mut_nt, generic_dna)
                sr = SeqRecord(mut_seq.translate(),
                               id='%s-%s-%s;freq=%f' % (lsp[3], lsp[1], lsp[4], freqs[i + 1]),
                               description='')
                all_muts.append(sr)

                '''
                for nt_pos in range(codon_pos * 3, (codon_pos + 1) * 3):
                    if pos - frame == nt_pos:
                        mut_cod += alt_nt
                        mut_cod += nt_ref[nt_pos + 1:]
                        break
                    else:
                        mut_cod += nt_ref[nt_pos]
                aa_ref_here = aa_ref[codon_pos]
                aa_alt_here = translation_table[mut_cod]
                
                if aa_ref_here != aa_alt_here:
                    allvars.append((aa_ref_here, codon_pos, aa_alt_here, freqs))
                    Dn += 1
                else:
                    Ds += 1
                '''
    SeqIO.write(all_muts, 'mutants.fasta', 'fasta')
#    print('Dn: %d\tDs: %d' % (Dn, Ds))
    return allvars


def main(ref_file='HXB2_pol_gene.fasta', bamfile='hq_2_cons_sorted.bam'):
    ''''''

    frame, nt_framed, aa_framed = find_frame(ref_file)

    # call variants with freebayes
    bamstem = '.'.join(bamfile.split('.')[:-1])
    cml = 'freebayes --min-alternate-total 5 --min-alternate-fraction 0.05 --ploidy 20'
    cml += ' --fasta-reference %s' % ref_file
    cml += ' %s > %s_fb.vcf' % (bamfile, bamstem)
    subprocess.call(cml, shell=True)

    filteredvars = parsevar(frame, nt_framed, aa_framed, '%s_fb.vcf' % bamstem) 
    return 'mutants.fasta'

if __name__ == '__main__':
    main()
