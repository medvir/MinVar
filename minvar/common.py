#!/usr/bin/env python3
"""Definitions, lists, resources, used throughout the tool."""

import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# from Table 1 of 10.1002/hep.20819, acc number of
# confirmed HCV genotypes/subtypes
acc_numbers = {
    '1a': ['M62321', 'M67463'],   # no gap in the polyprotein alignment, d=4.0%
    '1b': ['D90208', 'M58335'],   # no gap in the polyprotein alignment, d=5.3%
    '2a': ['D00944', 'AB047639'], # no gap in the polyprotein alignment, d=8.1%
    '2b': ['D10988', 'AB030907'], # no gap in the polyprotein alignment, d=4.2%
    '2c': ['D50409'],
    '2k': ['AB031663'],
    '3a': ['D17763', 'D28917'],  # no gap in the polyprotein alignment, d=5.5%
    '3b': ['D49374'],
    '3k': ['D63821'],
    '4a': ['Y11604'],
    # 4d
    '5a': ['Y13184'],
    '6a': ['Y12083'],
    'RF1_2k/1b': ['AY587845']
}

hcv_map = {}
for k, v in acc_numbers.items():
    for an in v:
        hcv_map[an] = k


# extract polyprotein
for s, v in hcv_map.items():
    cml = 'efetch -db nuccore -format gp -id %s -mode xml' % s
    cml += ' | xtract -pattern GBSeq -element GBSeq_locus -division GBFeature'
    cml += ' -block GBFeature_quals -section GBQualifier -if GBQualifier_name'
    cml += ' -equals translation -element GBQualifier_value'
    proc = subprocess.Popen(cml, shell=True, stdout=subprocess.PIPE,
                            universal_newlines=True)
    while True:
        line = proc.stdout.readline()
        if line:
            acn, sequence = line.split()
            if len(sequence) < 1000:
                break
            print('Writing %s, length: %d' % (acn, len(sequence)))
            sr = SeqRecord(Seq(sequence), id=acn, description='polyprotein')
            SeqIO.write(sr, '%s_%s.faa' % \
                        (v.replace('/', '_'), acn.split('.')[0]), 'fasta')
        else:
            break
org_dict = {
    'D50409.1': 'HCV',
    'D49374.1': 'HCV',
    'Y11604.1': 'HCV',
    'D63821.1': 'HCV',
    'D90208.1': 'HCV',
    'M58335.1': 'HCV',
    'M62321.1': 'HCV',
    'M67463.1': 'HCV',
    'D17763.1': 'HCV',
    'D28917.1': 'HCV',
    'D00944.1': 'HCV',
    'AB047639.1': 'HCV',
    'AB031663.1': 'HCV',
    'D10988.1': 'HCV',
    'AB030907.1': 'HCV',
    'AY587845.1': 'HCV',
    'CON_OF_CONS': 'HIV',
    'Mgroup': 'HIV',
    'CONSENSUS_A1': 'HIV',
    'A1.anc': 'HIV',
    'CONSENSUS_A2': 'HIV',
    'CONSENSUS_B': 'HIV',
    'B.anc': 'HIV',
    'CONSENSUS_C': 'HIV',
    'C.anc': 'HIV',
    'CONSENSUS_D': 'HIV',
    'CONSENSUS_F1': 'HIV',
    'CONSENSUS_F2': 'HIV',
    'CONSENSUS_G': 'HIV',
    'CONSENSUS_H': 'HIV',
    'CONSENSUS_01_AE': 'HIV',
    'CONSENSUS_02_AG': 'HIV',
    'CONSENSUS_03_AB': 'HIV',
    'CONSENSUS_04_CPX': 'HIV',
    'CONSENSUS_06_CPX': 'HIV',
    'CONSENSUS_08_BC': 'HIV',
    'CONSENSUS_10_CD': 'HIV',
    'CONSENSUS_11_CPX': 'HIV',
    'CONSENSUS_12_BF': 'HIV',
    'CONSENSUS_14_BG': 'HIV'
    }
