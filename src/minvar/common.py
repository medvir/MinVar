#!/usr/bin/env python3
"""Definitions, lists, resources, used throughout the tool."""

# from Table 1 of 10.1002/hep.20819 and
# from https://www.ncbi.nlm.nih.gov/projects/genotyping/view.cgi?db=3
# acc number of confirmed HCV genotypes/subtypes

acc_numbers = {
    '1a': ['M62321', 'M67463'],    # no gap in the polyprotein alignment, d=4.0%
    '1b': ['D90208', 'M58335'],    # no gap in the polyprotein alignment, d=5.3%
    '2a': ['D00944', 'AB047639'],  # no gap in the polyprotein alignment, d=8.1%
    '2b': ['D10988', 'AB030907'],  # no gap in the polyprotein alignment, d=4.2%
    '2c': ['D50409'],
    '2k': ['AB031663'],
    '3a': ['D17763', 'D28917'],  # no gap in the polyprotein alignment, d=5.5%
    '3b': ['D49374'],
    '3k': ['D63821'],
    '4a': ['Y11604'],
    '4d': ['DQ418786', 'FJ462437'],
    '5a': ['Y13184'],
    '6a': ['Y12083'],
    'RF1_2k/1b': ['AY587845']
}

hcv_map = {}
for k, v in acc_numbers.items():
    for an in v:
        hcv_map[an] = k

# dictionary with 1-based coordinates of genes on H77, 10.1002/hep.21377
h77_map = {}
for i in range(1, 192):
    h77_map[i] = 'C', i
for i in range(192, 384):
    h77_map[i] = 'E1', i - 191
for i in range(384, 747):
    h77_map[i] = 'E2', i - 383
for i in range(747, 810):
    h77_map[i] = 'p7', i - 746
for i in range(810, 1027):
    h77_map[i] = 'NS2', i - 809
for i in range(1027, 1658):
    h77_map[i] = 'NS3', i - 1028
for i in range(1658, 1712):
    h77_map[i] = 'NS4a', i - 1657
for i in range(1712, 1973):
    h77_map[i] = 'NS4b', i - 1711
for i in range(1973, 2421):
    h77_map[i] = 'NS5a', i - 1972
for i in range(2421, 3012):
    h77_map[i] = 'NS5b', i - 2420

# dictionary with 1-based coordinates of genes on consensus_B
consensus_B_map = {}
for i in range(1, 57):
    consensus_B_map[i] = 'GagPolTF', i
for i in range(57, 156):
    consensus_B_map[i] = 'protease', i - 56
for i in range(156, 596):
    consensus_B_map[i] = 'RT', i - 155
for i in range(596, 716):
    consensus_B_map[i] = 'RNase', i - 595
for i in range(716, 1004):
    consensus_B_map[i] = 'integrase', i - 715

org_dict = {
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
for k in hcv_map:
    org_dict[k] = 'HCV'
