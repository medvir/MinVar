#!/usr/bin/env python3
"""Parse the different output files and write a report in markdown, then convert it to pdf with markdown."""
import sys
import os
import csv
import logging
from pkg_resources import resource_filename

import pandas as pd

# from common import acc_numbers, hcv_map

# aminoacid one-letter code
aa_set = set('GPAVLIMCFYWHKRQNEDST')

# amminoacid sequences from files in db directory
# dn_dir = os.path.dirname(__file__)
# db_dir = os.path.abspath(os.path.join(dn_dir, 'db'))
# HCV_references = resource_filename(__name__, 'db/HCV/subtype_references.fasta')


# prot = \
#     list(SeqIO.parse(os.path.join(db_dir, 'protease.faa'), 'fasta'))[0]
# RT = list(SeqIO.parse(os.path.join(db_dir, 'RT.faa'), 'fasta'))[0]
# integrase = \
#     list(SeqIO.parse(os.path.join(db_dir, 'integrase.faa'), 'fasta'))[0]


def parse_drm():
    """Parse drug resistance mutations listed in files db/HIV/masterComments*.txt."""
    df_list = []
    genes = ['protease', 'RT', 'integrase']
    HIVdb_comment_files = [resource_filename(__name__, 'db/HIV/masterComments_%s.txt' % p)
                           for p in ('PI', 'RTI', 'INI')]
    for gene, drm_file in zip(genes, HIVdb_comment_files):
        # [('protease', 'masterComments_PI.txt'),
        #                        ('RT', 'masterComments_RTI.txt'),
        #                        ('integrase', 'masterComments_INI.txt')]:
        # drm_file = os.path.join(db_dir, 'HIV', drm_file_name)
        try:
            d1 = pd.read_table(drm_file, header=0,
                               names=['pos', 'mut', 'category', 'commented',
                                      'comment'])
        except pd.io.common.CParserError:  # integrase does not have commented column
            d1 = pd.read_table(drm_file, header=0,
                               names=['pos', 'mut', 'category', 'comment'])
        gs = pd.Series([gene] * len(d1))
        d1['gene'] = gs
        df_list.append(d1)
    df = pd.concat(df_list)
    return df


def parse_ras():
    """Parse position of RAS listed in file db/HCV/all_mutations_position.csv."""
    ras_file = resource_filename(__name__, 'db/HCV/all_mutations_position.csv')
    ras_positions = pd.read_csv(ras_file)
    ras_positions['CATEGORY'] = 'RAS'
    return ras_positions


def write_subtype_info(handle, subtype_file=None):
    """Write information on subtyping."""
    from operator import itemgetter
    md_header = 'Drug resistance mutations detected by NGS sequencing'
    md_header += '\n' + '=' * len(md_header) + '\n\n'
    md_header += 'Subtype inference with blast\n'
    md_header += '----------------------------\n'
    md_header += '|     subtype     | support [%] |\n'
    md_header += '|:{:-^15}:|:{:-^11}:|\n'.format('', '')
    if subtype_file:
        save_freq = {}
        with open(subtype_file) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for mtype, freq in spamreader:
                int_freq = int(round(100 * float(freq), 0))
                if int_freq >= 1:
                    save_freq[mtype] = save_freq.get(mtype, 0) + int_freq
        for k, v in sorted(save_freq.items(), key=itemgetter(1),
                           reverse=True):
            md_header += '|{: ^17}|{: ^13}|\n'.format(k, v)
    print(md_header, file=handle)


def write_header_HIV(handle, drms=None):
    """Write header to a file in markdown format."""

    md_header = """

Parsing mutations
-----------------

The list of annotated mutations was downloaded from HIVdb and includes:

"""
    for gene in ['protease', 'RT', 'integrase']:
        positions = set(drms[drms.gene == gene].pos.tolist())
        md_header += '- %d positions on %s\n' % (len(positions), gene)
    print(md_header, file=handle)


def write_header_HCV(handle, drms=None):
    """Write header to a file in markdown format."""

    md_header = """

Parsing mutations
-----------------

The list of known RAS includes:

"""
    for gene in ['NS3', 'NS5A', 'NS5B']:
        positions = set(drms[drms.gene == gene].pos.tolist())
        md_header += '- %d positions on %s\n' % (len(positions), gene)
    print(md_header, file=handle)


def aa_unpack(mut_string):
    """Helper function used to extract a set of amminoacids from the merged csv files."""
    if not mut_string.startswith('NOT'):
        return set(mut_string)
    return aa_set - set(mut_string.split()[1])


def parse_merged(mer_file):
    """Do this by hand because it was too complicated to achieve this functionality with panda alone."""
    with open(mer_file) as csvfile:
        reader = csv.DictReader(csvfile)
        mdf = pd.DataFrame(columns=reader.fieldnames)
        for row in reader:
            row['freq'] = float(row['freq'])
            row['pos'] = int(row['pos'])
            if row['mut_x'] in aa_unpack(row['mut_y']):
                mdf = mdf.append(row, ignore_index=True)
            elif row['mut_x'] == '-' and 'd' in aa_unpack(row['mut_y']):
                mdf = mdf.append(row, ignore_index=True)
            elif row['mut_x'] and row['mut_y'] == '':
                row['category'] = 'unannotated'
                mdf = mdf.append(row, ignore_index=True)
    return mdf


def write_contact_file(sample_id='unknown sample', version='unknown'):
    """Write tex file with contact information taken from ini file and sample_id in footer."""
    import configparser
    config = configparser.ConfigParser()
    config.read(os.path.expanduser('~/.minvar/contact.ini'))
    try:
        contact_dict = config['contact']
    except KeyError:
        contact_dict = {
            'unit': 'TO BE DEFINED',
            'phone': 'TO BE DEFINED',
            'fax': 'TO BE DEFINED',
            'email': 'TO BE DEFINED'
        }
    oh = open('contact.tex', 'w')
    oh.write(r'\fancyfoot[L]{sample: %s}' % sample_id + '\n')
    oh.write(r'\fancyfoot[C]{\texttt{built with MinVar version: %s}}' % version + '\n')
    oh.write(r'\begin{minipage}{0.5\textwidth}' + '\n')
    if 'logo' in contact_dict.keys():
        logo_file = os.path.expanduser(os.path.join('~/.minvar/', contact_dict['logo']))
        oh.write(r'\includegraphics[width=2in]{%s}' % logo_file + '\n')
    oh.write(r'\end{minipage}' + '\n')
    oh.write(r'\hfill' + '\n')
    oh.write(r'\begin{minipage}{0.45\textwidth}' + '\n')
    oh.write(r'\begin{tabular}{@{}r@{}}' + '\n')
    oh.write(r'\today \\[\normalbaselineskip]' + '\n')
    oh.write(r'%s\\' % contact_dict['unit'] + '\n')
    oh.write(r'phone: %s\\' % contact_dict['phone'] + '\n')
    oh.write(r'fax: %s\\' % contact_dict['fax'] + '\n')
    oh.write(r'email: \href{mailto:%s}{%s}\\' % (contact_dict['email'], contact_dict['email']) + '\n')
    oh.write(r'\end{tabular}' + '\n')
    oh.write(r'\end{minipage}' + '\n')

    r"""\fancyfoot[L]{John Smith \quad - \quad \href{mailto:john@smith.com}{{\small john@smith.com}}}

    \begin{minipage}{0.5\textwidth}
      \includegraphics[width=2in]{uzh_logo_e_pos}% Logo
    \end{minipage}
    \hfill
    \begin{minipage}{0.45\textwidth}
    \begin{tabular}{@{}r@{}}
      \today \\[\normalbaselineskip]
      Institute of Medical Virology \\
      Telephone + 41 44 63 42653  \\
      Fax + 41 44 63 44967 \\
      Email \href{mailto:med-virology@virology.uzh.ch}{med-virology@virology.uzh.ch}
    \end{tabular}
    \end{minipage}"""


def main(org=None, fastq=None, version='unknown', mut_file='final.csv', subtype_file='subtype_evidence.csv'):
    """What the main does."""
    import subprocess
    import shutil
    import re

    rh = open('report.md', 'w')
    write_subtype_info(rh, subtype_file)
    if org == 'HIV':
        resistance_mutations = parse_drm()
        write_header_HIV(rh, resistance_mutations)
    elif org == 'HCV':
        resistance_mutations = parse_ras()
        write_header_HCV(rh, resistance_mutations)

    logging.info('Reading mutations from %s', mut_file)
    mutation_detected = pd.read_csv(mut_file)
    logging.info('shape: %s', mutation_detected.shape)
    logging.info('Parsed DRM from database, shape: %s', str(resistance_mutations.shape))
    mpd = pd.merge(mutation_detected, resistance_mutations, how='left',
                   on=['gene', 'pos'])
    mpd.to_csv(path_or_buf='merged_muts_drm_annotated.csv', index=False)
    logging.info('Shape of raw merged is: %s', str(mpd.shape))

    if org == 'HIV':
        # too complicated with panda, do it by hand
        drms = parse_merged('merged_muts_drm_annotated.csv')
        logging.info('Shape of elaborated merged is: %s', str(drms.shape))
        # os.remove('merged_muts.csv')

        drms.drop(['commented', 'mut_y'], axis=1, inplace=True)
        drms.rename(columns={'mut_x': 'mut'}, inplace=True)

        cols = ['gene', 'pos', 'mut', 'freq', 'category']
        drms = drms[cols]
        drms.sort_values(by=['gene', 'pos', 'freq'], inplace=True,
                         ascending=[True, True, False])
        drms.to_csv('annotated_DRM.csv', index=False)

        for gene in ['protease', 'RT', 'integrase']:
            gene_muts = drms[drms.gene == gene]
            if gene_muts.shape[0] == 0:
                logging.info('No mutations on %s', gene)
                print('No mutations on %s' % gene, file=rh)
                print('-' * (len(gene) + 16), file=rh)
                print(file=rh)
                continue
            # sort was lost because
            gene_muts = gene_muts.sort_values(
                by=['gene', 'pos', 'freq'], ascending=[True, True, False])
            print('%s' % gene, file=rh)
            print('-'*len(gene), file=rh)
            h1 = '| position | mutation | frequency [%] |      category      |'
            print(h1, file=rh)
            h2 = '|:{:-^8}:|:{:-^8}:|:{:-^13}:|:{:-^18}:|'
            print(h2.format('', '', '', ''), file=rh)
            for index, row in gene_muts.iterrows():
                # # same pos mut tuple must give same annota, probably redundant
                # assert group['category'].nunique() == 1, \
                #     group['category'].describe()
                mut_cat = row['category']
                int_freq = int(round(100 * row['freq'], 0))
                print(
                    '|{: ^10}|{: ^10}|{: ^15}|{: ^20}|'.format(int(row['pos']),
                                                               row['mut'],
                                                               int_freq,
                                                               mut_cat),
                    file=rh)
            print('\n', file=rh)
    elif org == 'HCV':
        drms = pd.read_csv('merged_muts_drm_annotated.csv')
        drms = drms[drms['CATEGORY'] == 'RAS']
        # logging.info('Shape of annotated_mutations is: %s', str(mutation_detected.shape))

        if drms.shape[0] == 0:
            print('No mutations found', file=rh)
            print('-' * 26, file=rh)
            print('\n', file=rh)

        for gene in ['NS3', 'NS5A', 'NS5B']:
            gene_muts = drms[drms.gene == gene]
            if gene_muts.shape[0] == 0:
                logging.info('No RAS mutations on %s', gene)
                print('No RAS mutations on %s' % gene, file=rh)
                print('-' * (len(gene) + 16), file=rh)
                print(file=rh)
                continue

            print('%s' % gene, file=rh)
            print('-'*len(gene), file=rh)
            h1 = '| position | mutation | frequency [%] | type'
            print(h1, file=rh)
            h2 = '|:{:-^8}:|:{:-^8}:|:{:-^13}:|:{:-^10}:|'
            print(h2.format('', '', '', ''), file=rh)
            for index, row in gene_muts.iterrows():
                int_freq = int(round(100 * row['freq'], 0))
                category = row['CATEGORY']
                print(
                    '|{: ^10}|{: ^10}|{: ^15}|{: ^12}'.format(int(row['pos']), row['mut'], int_freq, category),
                    file=rh)
                del index
            print('\n', file=rh)
    rh.close()

    # copy template to current directory
    tmpl_file = resource_filename(__name__, 'db/template.tex')
    shutil.copy(tmpl_file, os.getcwd())
    # write contact.tex to be included in template
    fastq_base = os.path.basename(fastq)
    try:
        sample_id = re.search(r'(.*)_S\d*', fastq_base).group(1)
    except AttributeError:
        sample_id = 'unknown'
    write_contact_file(sample_id=sample_id, version=version)
    # convert to PDF with pandoc
    pand_cml = 'pandoc --template=./template.tex report.md -o report.pdf'
    logging.debug(pand_cml)
    subprocess.call(pand_cml, shell=True)
    os.remove('./template.tex')
    os.remove('contact.tex')


if __name__ == '__main__':
    main(org=sys.argv[1], fastq='xyz', version='unknown', mut_file='final.csv',
         subtype_file='subtype_evidence.csv')
