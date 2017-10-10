#!/usr/bin/env python3
'''Definitions, lists, resources used throughout the tool
'''
# from Table 1 of 10.1002/hep.20819, acc number of confirmed HCV genotypes/subtypes

acc_numbers = {
    '1a': ['M62321', 'M67463'],
    '1b': ['D90208', 'M58335'],
    '2a': ['D00944', 'AB047639'],
    '2b': ['D10988', 'AB030907'],
    '2c': ['D50409'],
    '2k': ['AB031663'],
    '3a': ['D17763', 'D28917'],
    '3b': ['D49374'],
    '3k': ['D63821'],
    '4a': ['Y11604'],
    '5a': ['Y13184'],
    '6a': ['Y12083'],
    'RF1_2k/1b': ['AY587845']
}

hcv_map = {}
for k, v in acc_numbers.items():
    for an in v:
        hcv_map[an] = k
