#!/usr/bin/env python3
"""Definitions, lists, resources, used throughout the tool."""

# from Table 1 of 10.1002/hep.20819 and
# from https://www.ncbi.nlm.nih.gov/projects/genotyping/view.cgi?db=3
# acc number of confirmed HCV genotypes/subtypes

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
h77_map_dict = {}
for i in range(1, 192):
    h77_map_dict[i] = 'C'
for i in range(192, 384):
    h77_map_dict[i] = 'E1'
for i in range(384, 747):
    h77_map_dict[i] = 'E2'
for i in range(747, 810):
    h77_map_dict[i] = 'p7'
for i in range(810, 1027):
    h77_map_dict[i] = 'NS2'
for i in range(1027, 1658):
    h77_map_dict[i] = 'NS3'
for i in range(1658, 1712):
    h77_map_dict[i] = 'NS4a'
for i in range(1712, 1973):
    h77_map_dict[i] = 'NS4b'
for i in range(1973, 2421):
    h77_map_dict[i] = 'NS5a'
for i in range(2421, 3012):
    h77_map_dict[i] = 'NS5b'


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
