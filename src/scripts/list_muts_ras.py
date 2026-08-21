#!/usr/bin/env python3
"""Parse RAS tables to extract specific dangerous mutations per genotype."""
import glob
import re
import pandas as pd

df_list = []
for f in glob.glob('NS*csv'):
    gene_name = f.split('_')[0]
    
    with open(f) as csvfile:
        muts = pd.read_csv(csvfile, delimiter=',')
        
        processed_rows = []
        for index, row in muts.iterrows():
            pos = int(row['aa_position'])
            new_row = {'gene': gene_name, 'pos': pos}
            
            # Keep the specific letters for each genotype column
            for col in muts.columns[1:]:
                val = str(row[col])
                dangerous_aas = set()
                
                if pd.notna(val) and val.strip() != '' and val != 'nan':
                    match = re.search(r'\d+([A-Z/]+)', val)
                    if match:
                        letters = match.group(1).split('/')
                        dangerous_aas.update(letters)
                
                if dangerous_aas:
                    new_row[col] = '/'.join(sorted(dangerous_aas))
                else:
                    new_row[col] = 'unknown'
                    
            processed_rows.append(new_row)
            
    df = pd.DataFrame(processed_rows)
    df_list.append(df)

mut_df = pd.concat(df_list)
mut_df = mut_df.sort_values(by=['gene', 'pos'])
mut_df.to_csv('all_mutations_position.csv', index=False)