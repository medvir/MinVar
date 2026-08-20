#!/usr/bin/env python3
"""Use this on RAS tables to have a list of all gene/positions listed."""
import glob
import re
import pandas as pd

df_list = []
for f in glob.glob('NS*csv'):
    gene_name = f.split('_')[0]
    pos_muts = {}
    
    with open(f) as csvfile:
        muts = pd.read_csv(csvfile, delimiter=',')
        
        for index, row in muts.iterrows():
            pos = int(row['aa_position'])
            dangerous_aas = set()
            
            # Iterate through all genotype columns (skip the first column 'aa_position')
            for col in muts.columns[1:]:
                val = str(row[col])
                # If the cell is not empty
                if pd.notna(val) and val.strip() != '' and val != 'nan':
                    # Regex: find the position number and capture the letters after it
                    # Example: from "L/F28M/V/S" it matches "28" and captures "M/V/S"
                    match = re.search(r'\d+([A-Z/]+)', val)
                    if match:
                        letters = match.group(1).split('/')
                        dangerous_aas.update(letters)
            
            if dangerous_aas:
                # Save as a slash-separated string like "C/D/E/G/H"
                pos_muts[pos] = '/'.join(sorted(dangerous_aas))

    df = pd.DataFrame({
        'pos': list(pos_muts.keys()),
        'dangerous_muts': list(pos_muts.values()),
        'gene': gene_name
    })
    df_list.append(df)

mut_df = pd.concat(df_list)
mut_df = mut_df.sort_values(by=['gene', 'pos'])
# Ensure columns are in the correct order
mut_df = mut_df[['gene', 'pos', 'dangerous_muts']]
mut_df.to_csv('all_mutations_position.csv', index=False)