#!/usr/bin/env python3.4

import sys
import os
import warnings
from pprint import pprint

import pandas as pd
from Bio import SeqIO

# aminoacid one-letter code
aa_set = set('GPAVLIMCFYWHKRQNEDST')

# amminoacid sequences from files in db directory
dn_dir = os.path.dirname(__file__)
db_dir = os.path.abspath(os.path.join(dn_dir, os.pardir, 'db'))
prot = \
    list(SeqIO.parse(os.path.join(db_dir, 'protease.faa'), 'fasta'))[0]
RT = list(SeqIO.parse(os.path.join(db_dir, 'RT.faa'), 'fasta'))[0]
integrase = \
    list(SeqIO.parse(os.path.join(db_dir, 'integrase.faa'), 'fasta'))[0]

def parse_drm():
    '''Parse drug resistance mutations listed in files db/*Variation.txt'''


    df_list = []
    for gene, drm_file_name in [('protease', 'masterComments_PI.txt'),
                                ('RT', 'masterComments_RTI.txt'),
                                ('integrase', 'masterComments_INI.txt')]:
        drm_file = os.path.join(db_dir, drm_file_name)
        try:
            d1 = pd.read_table(drm_file, header=0, names=['pos', 'mut', 'category',
                               'commented', 'comment'])
        except:  # integrase does not have commented column
            d1 = pd.read_table(drm_file, header=0, names=['pos', 'mut', 'category',
                               'comment'])
        gs = pd.Series([gene] * len(d1))
        d1['gene'] = gs
        df_list.append(d1)
    df = pd.concat(df_list)
    return df


def parse_region(hap1):
    '''Whether the region comes from protease, RT or integrase'''

    import Alignment
    matched_res = []

    for gs in [prot, RT, integrase]:
        alout = '%s.tmp' % gs.id
        Alignment.needle_align('asis:%s' % hap1, 'asis:%s' % str(gs.seq),
                               alout)
        ald = Alignment.alignfile2dict([alout])
        os.remove(alout)
        ali = ald['asis']['asis']
        ali.summary()
        ratio1 = float(ali.ident) / (ali.ident + ali.mismatches)
        if ratio1 > 0.75:
            matched_res.append(gs.id)

    return matched_res


def write_header(handle, matched_types=None):
    '''Write header to a file in markdown format'''
    md_header = 'Drug resistance mutations detected by NGS sequencing'
    mc = len(md_header)
    md_header += '\n' + '=' * len(md_header) + '\n\n'
    md_header +='Subtype inference with blast\n'
    md_header +='----------------------------\n'
    md_header += '|    subtype    | support [%] |\n'
    md_header += '|:{:-^13}:|:{:-^11}:|\n'.format('', '', '')
    if matched_types:
        for mtype, freq in matched_types:
            int_freq = int(round(100 * freq, 0))
            md_header += '|{: ^15}|{: ^13}|\n'.format(mtype, int_freq)
    md_header += '''

Parsing mutations
-----------------

The list of mutations was downloaded from HIVdb and includes:

- xyz positions on protease
- zyx positions on RT
- abc positions on integrase.
'''
    print(md_header, file=handle)


def parse_com_line():
    '''argparse, since optparse deprecated starting from python 2.7'''
    import argparse

    parser = argparse.ArgumentParser(description='parse DRM mutations')
    parser.add_argument('-a', '--haplotypes', dest='haps',
                        help='file with amino acid haplotypes')

    args = parser.parse_args()
    return args


def aa_unpack(mut_string):
    if not mut_string.startswith('NOT'):
        return set(mut_string)
    else:
        return aa_set - set(mut_string.split()[1])


def parse_merged(mer_file):
    '''This is done by hand because it was too complicated to achieve
    this functionality with panda alone'''
    import csv
    mdf = []
    with open(mer_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            row['freq'] = float(row['freq'])
            row['pos'] = int(row['pos'])
            if row['mut_x'] in aa_unpack(row['mut_y']):
                mdf.append(row)
            elif row['mut_x'] and row['mut_y'] == '':
                row['category'] = 'unannotated'
                mdf.append(row)
    return pd.DataFrame(mdf)

def main(hap_file, matched_types=None):
    '''What does the main do?'''
    import subprocess
    import pickle

    # Appending to path in order to find Alignment
    sys.path.append(dn_dir)

    rh = open('report.md', 'w')
    write_header(rh, matched_types)

    if not hap_file:
        args = parse_com_line()
        hap_file = args.haps

    haps = list(SeqIO.parse(hap_file, 'fasta'))
    print('Searching genes for match. This only reports genes covered at least \
at 75%.\nIndividual aminoacid calls might be present. \n', file=rh)
    matched_region = parse_region(str(haps[0].seq))
    if not matched_region:
        print('None found: mutations not in protease/RT/integrase.\n', file=rh)
        sys.exit('Haplotypes not in protease/RT/integrase')

    for g in ['protease', 'RT', 'integrase']:
        mt = 'yes' if g in matched_region else 'no'
        print('**%s: %s**  ' % (g, mt), file=rh)
    print ('', file=rh)

    stem = '.'.join(hap_file.split('.')[:-1])
    csv_file = 'mutations.csv'

    print('Parsing DRM from database', file=sys.stderr)
    resistance_mutations = parse_drm()
    print('Shape is: ', resistance_mutations.shape)

    print('Reading mutations from %s' % csv_file, file=sys.stderr)
    mutation_detected = pd.read_csv(csv_file)
    print('Shape is: ', mutation_detected.shape)

    mpd = pd.merge(mutation_detected, resistance_mutations, how='left',
                   on=['gene', 'pos'])
    mpd.to_csv(path_or_buf='merged_muts.csv')
    print('Shape of raw merged is: ', mpd.shape)

    # too complicated with panda, do it by hand
    drms = parse_merged('merged_muts.csv')
    print('Shape of merged is: ', drms.shape)
    os.remove('merged_muts.csv')

    drms.drop(['', 'commented', 'mut_y'], axis=1, inplace=True)
    drms.rename(columns={'mut_x': 'mut'}, inplace=True)

    cols = ['gene', 'pos', 'mut', 'freq', 'category']
    drms = drms[cols]
    drms.sort(cols[:3], inplace=True)
    drms.to_csv('annotated_mutations.csv', index=False)

    for gene in ['protease', 'RT', 'integrase']:
        gene_muts = drms[drms.gene == gene]
        if gene_muts.shape[0] == 0:
            print('No mutations on ', gene, file=sys.stderr)
            continue

        grouped = gene_muts.groupby(['pos', 'mut'])
        print('%s' % gene, file=rh)
        print('-'*len(gene), file=rh)
        print('| position | mutation | frequency [%] |      category      |',
              file=rh)
        print('|:{:-^8}:|:{:-^8}:|:{:-^13}:|:{:-^18}:|'.format('', '', '', ''),
              file=rh)
        for name, group in grouped:
            # same pos mut tuple must give same annotation, probably redundant
            assert group['category'].nunique() == 1, group['category'].describe()
            mut_cat = group['category'].unique()[0]
            int_freq = int(round(100 * group['freq'].sum(), 0))
            print('|{: ^10}|{: ^10}|{: ^15}|{: ^20}|'.format(int(name[0]),
                                                     name[1],
                                                     int_freq,
                                                     mut_cat), file=rh)
        print('\n', file=rh)
    rh.close()
    # convert to PDF with pandoc
    tmpl_file = os.path.abspath(os.path.join(dn_dir, 'template.tex'))
    pand_cml = 'pandoc --template={} report.md -o report.pdf'.format(tmpl_file)
    print(pand_cml, file=sys.stderr)
    subprocess.call(pand_cml, shell=True)


if __name__ == '__main__':
    args = parse_com_line()
    main(args.haps)
