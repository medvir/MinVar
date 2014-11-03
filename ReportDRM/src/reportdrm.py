#!/usr/bin/env python3.4
'''Doctext to add
'''

import sys
import os
import warnings

import pandas as pd
from Bio import SeqIO

# amminoacid sequences from files in db directory
dn_dir = os.path.dirname(__file__)
db_dir = os.path.abspath(os.path.join(dn_dir, os.pardir, 'db'))
prot = \
    list(SeqIO.parse(os.path.join(db_dir, 'protease.fasta'), 'fasta'))[0]
RT = list(SeqIO.parse(os.path.join(db_dir, 'RT.fasta'), 'fasta'))[0]
integrase = \
    list(SeqIO.parse(os.path.join(db_dir, 'integrase.fasta'), 'fasta'))[0]

def parse_drm():
    '''Parse drug resistance mutations listed in files db/*Variation.txt'''


    df = pd.DataFrame(columns=['pos', 'mut', 'gene'])
    for gene, drm_file_name in [('protease', 'PRVariation.txt'),
                                ('RT', 'RTVariation.txt'),
                                ('integrase', 'INVariation.txt')]:
        drm_file = os.path.join(db_dir, drm_file_name)
        for l in open(drm_file):
            if not l.startswith('#'):
                pos, mut = l.strip().split('|')
                for m in mut:
                    df = df.append({'gene': gene, 'pos': int(pos), 'mut': m},
                                   ignore_index=True)
    return df


def parse_mutations(haplos):
    '''returns a list of mutations with frequencies'''
    import Alignment

    df = pd.DataFrame(columns=['haplotype', 'wt', 'pos', 'mut', 'freq', 'gene'])

    freqs = []
    dict_list = []
    for h in haplos:
        freq = float(h.description.split('=')[1])
        for mreg in [prot, RT, integrase]:
            Alignment.needle_align('asis:%s' % str(mreg.seq),
                                   'asis:%s' % str(h.seq), 'h.tmp')
            alhr = Alignment.alignfile2dict(['h.tmp'])

            os.remove('h.tmp')
            alih = alhr['asis']['asis']
            alih.summary()
            if 3 * alih.mismatches > alih.ident:
                continue
            wt, mut = alih.seq_a, alih.seq_b

            i = 0
            end_wt, end_mut = len(wt.rstrip('-')), len(mut.rstrip('-'))
            end_pos = min(end_wt, end_mut)
            for z in list(zip(wt, mut))[:end_pos]:
                i += z[0] != '-'
                if z[0] != z[1] and i:
                    mutation = '%s%d%s' % (z[0], i, z[1])
                    dict_here = {'haplotype': h.id, 'gene': mreg.id, 'wt': z[0],
                                 'pos': i, 'mut': z[1], 'freq': freq}

                    dict_list.append(dict_here)

    df = df.append(dict_list)
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


def write_header(handle, match_id=None):
    '''Write header to a file in markdown format'''
    md_header = 'Drug resistance mutations detected by NGS sequencing'
    mc = len(md_header)
    md_header += '\n' + '=' * len(md_header) + '\n\n'
    if match_id:
        md_header += 'Subtype detected: {}'.format(match_id)
    md_header += '''

The list of mutations was downloaded from HIVdb and includes:

- 68 positions on protease
- 344 positions on RT
- 134 positions on integrase.
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

def main(hap_file, match_id=None):
    '''What does the main do?'''
    import subprocess
    import pickle

    # Appending to path in order to find Alignment
    sys.path.append(dn_dir)
    

    rh = open('report.md', 'w')
    write_header(rh)

    if not hap_file:
        args = parse_com_line()
        hap_file = args.haps

    haps = list(SeqIO.parse(hap_file, 'fasta'))
    print('Searching genes for match. This only reports genes covered at least \
at 75 %.\nIndividual aminoacid calls might be present. \n', file=rh)
    matched_region = parse_region(str(haps[0].seq))
    if not matched_region:
        print('None found: mutations not in protease/RT/integrase.\n', file=rh)
        sys.exit('Haplotypes not in protease/RT/integrase')

    for g in ['protease', 'RT', 'integrase']:
        mt = 'yes' if g in matched_region else 'no'
        print('**%s: %s**  ' % (g, mt), file=rh)
    print ('', file=rh)

    stem = '.'.join(hap_file.split('.')[:-1])
    json_file = stem + '.json'

    if not os.path.exists(json_file):
        print('Aligning %s' % hap_file, file=sys.stderr)

        resistance_mutations = parse_drm()
        mutation_detected = parse_mutations(haps)

        mpd = pd.merge(mutation_detected, resistance_mutations,
                       on=['gene', 'pos', 'mut'])
        # drop haplotype and wild type
        mpd = mpd.iloc[:, 2:]
        mpd.to_json(path_or_buf=json_file)
    else:
        print('JSON file', json_file, 'found', file=sys.stderr)
        mpd = pd.read_json(json_file)

    for gene in ['protease', 'RT', 'integrase']:
        gene_muts = mpd[mpd.gene == gene]
        if gene_muts.shape[0] == 0:
            continue
        grouped = gene_muts.groupby(['pos', 'mut'])

        print('%s' % gene, file=rh)
        print('-'*len(gene), file=rh)
        print('| position | mutation | frequency [%] |', file=rh)
        print('|:{:-^8}:|:{:-^8}:|:{:-^13}:|'.format('', '', ''), file=rh)
        for name, group in grouped:
            int_freq = int(round(100 * group['freq'].sum(), 0))
            print('|{: ^10}|{: ^10}|{: ^15}|'.format(int(name[0]),
                                                     name[1],
                                                     int_freq),
                  file=rh)
        print('\n', file=rh)

    # convert to PDF with pandoc
    tmpl_file = os.path.abspath(os.path.join(dn_dir, os.pardir, 'src/template.tex'))
    subprocess.call('pandoc --template={} report.md -o report.pdf'.format(tmpl_file),
                    shell=True)


if __name__ == '__main__':
    args = parse_com_line()
    main(args.haps)
