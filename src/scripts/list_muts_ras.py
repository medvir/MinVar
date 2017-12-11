#!/usr/bin/env python3
"""Use this on RAS tables to have a list of all gene/positions listed."""
import glob
import pandas as pd

df_list = []
for f in glob.glob('NS*csv'):
    gene_name = f.split('_')[0]
    with open(f) as csvfile:
        muts = pd.read_csv(csvfile, delimiter=',')
        df = pd.DataFrame({'aa_position': list(set(muts.aa_position)), 'gene': gene_name})
        df_list.append(df)
mut_df = pd.concat(df_list)
mut_df.rename(columns={'aa_position': 'pos'}, inplace=True)
mut_df = mut_df.sort_values(by=['gene', 'pos'])
mut_df = mut_df[['gene', 'pos']]
mut_df.to_csv('all_mutations_position.csv', index=False)
