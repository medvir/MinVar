#!/usr/bin/env python
from pprint import pprint
import subprocess

'''' Download ref sequences from ICTV maintained file
download with
wget -O HCV_ref.fasta https://talk.ictvonline.org/ictv_wikis/flaviviridae/hepacivirus/m/hepacivirus-files/6789/download
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
    '4a': ['Y11604']
}

genotype = {}

acc_list = []
for k, v in acc_numbers.items():
    for n in v:
        acc_list.append(n)
        genotype[n] = k
pprint(genotype)

cml = 'efetch -db nuccore -format fasta -id \"%s\" > HCV_ref.fasta' % ','.join(acc_list)
#subprocess.call(cml, shell=True)
