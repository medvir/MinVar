#!/usr/bin/env python
'''' Download ref sequences from ICTV maintained file with
wget -O HCV_ref.fasta https://talk.ictvonline.org/ictv_wikis/flaviviridae/hepacivirus/m/hepacivirus-files/6789/download
or interrogate NCBI with efetch
'''

import os
import sys
import subprocess
from pprint import pprint

# manipulate path to import functions
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(1, parent_dir)
mod = __import__('minvar')
sys.modules["minvar"] = mod
from minvar.common import acc_numbers

genotype = {}
acc_list = []
for k, v in acc_numbers.items():
    for n in v:
        acc_list.append(n)
        genotype[n] = k
pprint(genotype)

cml = 'efetch -db nuccore -format fasta -id \"%s\" > HCV_ref.fasta' % ','.join(acc_list)
subprocess.call(cml, shell=True)
