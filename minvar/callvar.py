#!/usr/bin/env python3
'''Calls lofreq to produce a vcf file'''

import sys
import subprocess
import os
import warnings
import pandas as pd

import logging

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

# manipulate path to import functions
dn_dir = os.path.dirname(os.path.abspath(__file__))
os.sys.path.insert(1, dn_dir)
import Alignment

RAW_DEPTH_THRESHOLD = 50
MIN_FRACTION = 0.015
MAPPING_QUALITY_THRESHOLD = 20

try:
    CPUS = max(1, os.cpu_count())
except TypeError:
    CPUS = 1

# amminoacid sequences from files in db directory
db_dir = os.path.abspath(os.path.join(dn_dir, 'db'))
B_pol_nt = \
    list(SeqIO.parse(os.path.join(db_dir, 'consensus_B.fna'), 'fasta'))[0]
B_pol_aa = \
    list(SeqIO.parse(os.path.join(db_dir, 'consensus_B.faa'), 'fasta'))[0]

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


#def recalibrate_qualities(ref_file, bamfile, platform="454"):
def recalibrate_qualities(ref_file, bamfile, platform="ILLUMINA"):
    '''Invoke GATK BaseRecalibrator, also calling some high confidence variants.
    Follows vipr by Andreas Wilm'''

    bamstem, bamext = os.path.splitext(bamfile)
    assert bamext == '.bam', bamext
    recalfile = '%s_recal.bam' % bamstem
    if os.path.exists(recalfile):
        logging.debug('file %s exists, not overwriting it' % recalfile)
        return recalfile
    # count mapped reads
    cml = "samtools view -c -F 4 %s" % bamfile
    mapped_reads = int(subprocess.check_output(cml, shell=True))

    # sample 10k reads and run lofreq to detect high confidence variants
    fraction_needed = min(1.0, round(float(10000) / mapped_reads, 3))
    cml = 'samtools view -s %f -F 4 -h -b -o subsampled.bam %s' % (fraction_needed, bamfile)
    subprocess.call(cml, shell=True)
    cml = 'samtools index subsampled.bam'
    subprocess.call(cml, shell=True)

    cml = 'lofreq call-parallel --pp-threads %d -f %s -o tmp.vcf subsampled.bam' % (min(6, CPUS), ref_file)
    subprocess.call(cml, shell=True)
    cml = 'lofreq filter -v 200 -V 2000 -a 0.40 -i tmp.vcf -o known.vcf'
    subprocess.call(cml, shell=True)
    os.remove('tmp.vcf')

    # need to add group information to reads
    print("@RG\tID:minvar\tSM:haploid\tLB:ga\tPL:%s" % platform, file=open('rg.txt', 'w'))
    cml = "samtools view -h %s |  cat rg.txt - | " % bamfile
    cml += "awk '{ if (substr($1,1,1)==\"@\") print; else printf \"%s\\tRG:Z:minvar\\n\",$0; }' | "
    cml += "samtools view -u - > grouped.bam"
    subprocess.call(cml, shell=True)
    os.remove('rg.txt')
    cml = 'samtools index grouped.bam'
    subprocess.call(cml, shell=True)

    refstem = os.path.splitext(ref_file)[0]
    cml = 'java -jar /usr/local/picard-tools/picard.jar CreateSequenceDictionary R=%s O=%s.dict' % (ref_file, refstem)
    subprocess.call(cml, shell=True)

    # first pass to parse features
    cml = 'java -jar /usr/local/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator --maximum_cycle_value 600'
    cml += ' -I grouped.bam -l INFO -R %s -o recal_data.grp -knownSites known.vcf' % ref_file
    subprocess.call(cml, shell=True)

    # second pass to recalibrate
    cml = 'java -jar /usr/local/GATK/GenomeAnalysisTK.jar -T PrintReads'
    cml += ' -R %s -I grouped.bam -BQSR recal_data.grp -o %s' % (ref_file, recalfile)
    subprocess.call(cml, shell=True)

    return recalfile


def indelqual(ref_file, bamfile):
    '''Use lofreq indelqual to insert indel qualities into bamfile; this allows
    calling indels and is an alternative to GATK BaseRecalibrator
    '''

    bamstem, bamext = os.path.splitext(bamfile)
    assert bamext == '.bam', bamext
    recalfile = '%s_recal.bam' % bamstem
    cml = 'lofreq indelqual --dindel -f %s -o %s %s' % (ref_file, recalfile, bamfile)
    subprocess.call(cml, shell=True)
    # needs reindexing
    subprocess.call('samtools index %s' % recalfile, shell=True)

    return recalfile


def call_variants(ref_file=None, bamfile=None, parallel=True, n_regions=8, caller='lofreq'):
    '''Wrapper to call lofreq (also other tools originally)'''
    bamstem, bamext = os.path.splitext(bamfile)  # '.'.join(bamfile.split('.')[:-1])
    assert bamext == '.bam', bamfile + bamext

    # index reference
    cml = 'samtools faidx %s' % ref_file
    subprocess.call(cml, shell=True)

    # call minority variants with freebayes, lofreq or samtools/bcftools

    # freebayes here
    if parallel and caller == 'freebayes':
        # first compute regions
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
        cml += ' --min-alternate-count 10 --min-coverage %s' % RAW_DEPTH_THRESHOLD
        cml += ' --min-alternate-fraction %f --ploidy 20' % MIN_FRACTION
        cml += ' --haplotype-length 50'
        cml += ' --fasta-reference %s %s > %s.vcf' % (ref_file, bamfile, bamstem)
        subprocess.call(cml, shell=True)
        return '%s.vcf' % bamstem

    elif not parallel and caller == 'freebayes':
        cml = 'freebayes'
        cml += ' --min-alternate-count 10 --min-coverage %s' % RAW_DEPTH_THRESHOLD
        cml += ' --min-alternate-fraction %f --ploidy 20 --pooled-continuous' % MIN_FRACTION
        cml += ' --haplotype-length 50'
        cml += ' --fasta-reference %s %s > %s.vcf' % (ref_file, bamfile, bamstem)
        subprocess.call(cml, shell=True)
        return '%s.vcf' % bamstem

    # lofreq here
    elif parallel and caller == 'lofreq':
        cml = 'lofreq call-parallel --pp-threads %d --call-indels -f %s -o tmp.vcf %s' % \
            (min(6, CPUS), ref_file, bamfile)
        subprocess.call(cml, shell=True)
    elif not parallel and caller == 'lofreq':
        cml = 'lofreq call -f %s -o tmp.vcf %s' % (ref_file, bamfile)
        subprocess.call(cml, shell=True)

    # samtools here
    elif caller == 'samtools':
        # write a file with read group and sample info
        print("@RG\tID:minvar\tSM:haploid\tLB:ga", file=open('rg.txt', 'w'))
        cml = "samtools view -h %s" % bamfile
        cml += " | cat rg.txt - | awk \'{ if (substr($1,1,1)==\"@\") print; else printf \"%s\\tRG:Z:minvar\\n\",$0; }\'"
        cml += " | samtools view -u - > gr.bam"
        subprocess.call(cml, shell=True)
        # now mpileup and then bcftools with samples info
        cml = 'samtools mpileup -B -Q 20 -uf %s gr.bam > gr.mpileup' % ref_file
        subprocess.call(cml, shell=True)
        print("haploid\t1", file=open('samples.txt', 'w'))
        cml = 'bcftools call -S samples.txt -m -v gr.mpileup >  %s.vcf' % bamstem
        subprocess.call(cml, shell=True)
        return '%s.vcf' % bamstem

    # lofreq still needs filtering
    if caller == 'lofreq':
        cml = 'lofreq filter -i tmp.vcf -o %s.vcf -v %d -a %f' % \
            (bamstem, RAW_DEPTH_THRESHOLD, MIN_FRACTION)
        subprocess.call(cml, shell=True)
        return '%s.vcf' % bamstem


def main(ref_file=None, bamfile=None, parallel=True, n_regions=6,
         caller='lofreq', recalibrate=True):
    '''What does a main do?'''

    # call variants
    bamstem, bamext = os.path.splitext(bamfile)
    assert bamext == '.bam', bamext

    called_file = '%s.vcf' % bamstem
    if os.path.exists(called_file):
        logging.info('vcf file exists, not calling new variants')
        return called_file

    if recalibrate:
        logging.info('recalibrating base qualities')
        recalfile = recalibrate_qualities(ref_file, bamfile)
    else:
        logging.info('base qualities will not be recalibrated')
        if(caller == 'lofreq'):
            logging.info('indel qualities introduced with lofreq')
            recalfile = indelqual(ref_file, bamfile)
        else:
            logging.info('indel qualities will not be introduced')
            recalfile = bamfile

    logging.info('calling variants with %s' % caller)
    called_file = call_variants(ref_file, recalfile, parallel, n_regions, caller)
    called_bam = recalfile

    return called_file, called_bam

if __name__ == '__main__':
    main()
