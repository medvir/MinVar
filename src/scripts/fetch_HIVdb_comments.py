#!/usr/bin/env python3
"""Access DR Comments from HIVdb and arrange them in tsv files.
When masterComments files are updated, mastercomments_version in src/minvar/common.py must be updated too.
"""
import pandas as pd
import re

classes = ['PI', 'NRTI', 'NNRTI', 'INSTI']  # they are split over NRTI/NNRTI
td = {}
for c in classes:
    table = pd.read_html(r'https://hivdb.stanford.edu/dr-summary/comments/%s/' % c)[0]
    table.rename(columns={'Comment/Mutation Type': 'Category'}, inplace=True)
    # extract position/amminoacid
    table['POSITION'], table['AA'] = zip(
        *table.apply(lambda row: re.search('(\d*)(\w*)', row['Condition']).group(1, 2), axis=1))
    td[c] = table

# merge RTI
td['RTI'] = pd.concat([td['NRTI'], td['NNRTI']])
del td['NRTI']
del td['NNRTI']
# rename INSTI to INI
td['INI'] = td.pop('INSTI')

for k, v in td.items():
    v.to_csv('masterComments_%s.txt' % k, index=False, sep='\t', columns=['POSITION', 'AA', 'Category', 'Comment'])
