#!/usr/bin/env python3
"""Parse the different output files and write a report in markdown, then convert it to pdf with markdown."""
import sys
import os
import configparser
import csv
import logging
import shlex
import subprocess
from pkg_resources import resource_filename
from Bio import SeqIO

import pandas as pd

dn_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if __name__ == '__main__':
    if __package__ is None:
        os.sys.path.insert(1, dn_dir)
        mod = __import__('minvar')
        sys.modules["minvar"] = mod
        from common import MIN_FRACTION, RAW_DEPTH_THRESHOLD, drug_names, mastercomments_version
else:
    from .common import MIN_FRACTION, RAW_DEPTH_THRESHOLD, drug_names, mastercomments_version

# aminoacid one-letter code
aa_set = set('GPAVLIMCFYWHKRQNEDST')

# mutations below this threshold do not contribute to drug prediction via HIVdb
runpars = configparser.ConfigParser()
runpars.read(os.path.expanduser('~/.minvar/runpars.ini'))
try:
    sierra_threshold = float(runpars['sierra']['threshold'])
except KeyError:
    sierra_threshold = 0.2

# colour cells according to susceptibility
cell_colour = {
    'Susceptible': 'resistance2',
    'Potential Low-Level Resistance': 'resistance3',
    'Low-Level Resistance': 'resistance4',
    'Intermediate Resistance': 'resistance5',
    'High-Level Resistance': 'resistance6'
}

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
        d1 = pd.read_table(drm_file, header=0, names=['pos', 'mut', 'category', 'comment'])
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
    try:
    #if os.path.exists(subtype_file):
        md_header = 'Drug resistance mutations detected by NGS sequencing'
        md_header += '\n' + '=' * len(md_header) + '\n\n'
        save_freq = {}
        with open(subtype_file) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for mtype, freq in spamreader:
                int_freq = int(round(100 * float(freq), 0))
                if int_freq >= 1:
                    save_freq[mtype] = save_freq.get(mtype, 0) + int_freq
        top_gt, support = sorted(save_freq.items(), key=itemgetter(1), reverse=True)[0]
        if support >= 50:
            md_header += 'Inferred subtype: %s (blast support: %d %%)\n' % (top_gt, support)
            md_header += '-----------------------------------------\n'
        else:
            md_header += 'Subtype not inferred: (blast support too low)\n'
            md_header += '---------------------------------------------\n'

    except FileNotFoundError:
        md_header = '\n'
    print(md_header, file=handle)

def write_sierra_results(handle, mut_file):
    """Run sierrapy patterns results."""

    import json
    mutations = pd.read_csv(mut_file)
    mutations = mutations[mutations['freq'] >= sierra_threshold]
    mutations = mutations[mutations['gene'] != 'GagPolTF']
    mutations = mutations[mutations['gene'] != 'RNase']
    mutations = mutations[mutations['wt'] != '*']
    gmap = {'protease': 'PR', 'RT': 'RT', 'integrase': 'IN'}
    gmap_inverse = {'PR': 'protease', 'IN': 'integrase'}
    mutations['pattern'] = mutations.apply(
        lambda m: '%s:%s%d%s' % (gmap[m['gene']], m['wt'], m['pos'], m['mut']), axis=1)
    ptn = ' + '.join(mutations['pattern'])
    with open('pattern.txt', 'w') as h:
        h.write(ptn + '\n')
    cml = shlex.split('sierrapy patterns pattern.txt -o o.json')
    logging.debug(cml)
    subprocess.call(cml)
    with open('o.json') as h:
        patterns = json.load(h)
    os.remove('o.json')
    os.remove('pattern.txt')
    assert len(patterns) == 1
    pattern = patterns[0]
    version = pattern['drugResistance'][0]['version']['text']
    pubdate = pattern['drugResistance'][0]['version']['publishDate']
    print('Drug Resistance Interpretation', file=handle)
    print('==============================\n', file=handle)
    print('Stanford HIVdb version %s, pubdate %s.' % (version, pubdate), file=handle)
    print('Mutations below %d%% were not included.\n' % (100 * sierra_threshold), file=handle)

    for dr in pattern['drugResistance']:
        gene_name = dr['gene']['name']
        comments_set = set()
        print(gmap_inverse.get(gene_name, gene_name), file=handle)
        print('-' * len(gene_name) + '\n', file=handle)
        print('| class |         name         | score|      assessment      |               mutations                |',
              file=handle)
        print('|:{:-^5}:|:{:-^20}:| {:-^4}:|:{:-^20}:|:{:-^39}|'.format('', '', '', '', ''), file=handle)

        for drugscore in dr['drugScores']:
            drugClass = drugscore['drugClass']['name']
            drug = drug_names.get(drugscore['drug']['name'], drugscore['drug']['name'])
            drug += ' (%s)' % drugscore['drug']['displayAbbr']
            all_muts = ' + '.join((partial['mutations'][0]['text'] for partial in drugscore['partialScores']))

            colour = cell_colour.get(drugscore['text'], 'white')
            print('|{: ^7}|{: ^22}|{: ^6}|\\cellcolor{{{}}}{: ^22}|{: ^40}|'.format(
                drugClass, drug, drugscore['score'], colour, drugscore['text'], all_muts), file=handle)
            for partial in drugscore['partialScores']:
                for mutations in partial['mutations']:
                    for pmm in mutations['comments']:
                        comments_set.add(pmm['text'])
        print('', file=handle)
        for cm in comments_set:
            print('- %s' % cm, file=handle)

        print('', file=handle)
            # for partial in drugscore['partialScores']:
            #  if len(partial['mutations']) > 1:
            #      pprint(partial)
            # assert len(partial['mutations']) == 1, partial['mutations']
    print('\\newpage', file=handle)


def write_ambig_score(handle):
    """Read ambiguous sequence, compute and write ambiguity score."""

    ambi = list(SeqIO.parse('cns_ambiguous.fasta', 'fasta'))[0].seq
    ambi = str(ambi).replace('N', '')
    ambi_score = float(sum((1 for nt in ambi if not nt in set(['A', 'C', 'G', 'T']))))
    ambi_score /= len(ambi)
    ambi_score *= 100
    print('\nAmbiguity score: %3.1f %%' % ambi_score, file=handle)
    print('------------------------\n', file=handle)

def write_header_HIV(handle, drms=None):
    """Write header to a file in markdown format."""

    md_header = """\
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
    oh.write(r'\fancyfoot[C]{\texttt{MinVar version: %s}}' % version + '\n')
    oh.write(r'\begin{minipage}{0.5\textwidth}' + '\n')
    if 'logo' in contact_dict.keys():
        logo_file = os.path.expanduser(os.path.join('~/.minvar/', contact_dict['logo']))
        oh.write(r'\includegraphics[width=2in]{%s}' % logo_file + '\n')
    oh.write(r'\end{minipage}' + '\n')
    oh.write(r'\hfill' + '\n')
    oh.write(r'\begin{minipage}{0.45\textwidth}' + '\n')
    oh.write(r'\begin{tabular}{@{}r@{}}' + '\n')
    #oh.write(r'\today \\[\normalbaselineskip]' + '\n')
    oh.write(r'%s\\' % contact_dict['unit'] + '\n')
    oh.write(r'phone: %s\\' % contact_dict['phone'] + '\n')
    oh.write(r'fax: %s\\' % contact_dict['fax'] + '\n')
    oh.write(r'email: \href{mailto:%s}{%s}\\' % (contact_dict['email'], contact_dict['email']) + '\n')
    oh.write(r'\end{tabular}' + '\n')
    oh.write(r'\end{minipage}' + '\n')
    oh.write(r'\vspace{1cm}' + '\n')

    r"""\fancyfoot[L]{John Smith \quad - \quad \href{mailto:john@smith.com}{{\small john@smith.com}}}

    \begin{minipage}{0.5\textwidth}
      \includegraphics[width=2in]{uzh_logo_e_pos}% Logo
    \end{minipage}
    \hfill
    \begin{minipage}{0.45\textwidth}
    \begin{tabular}{@{}r@{}}
      %\DTMnow \\[\normalbaselineskip]
      Institute of Medical Virology \\
      Telephone + 41 44 63 42653  \\
      Fax + 41 44 63 44967 \\
      Email \href{mailto:med-virology@virology.uzh.ch}{med-virology@virology.uzh.ch}
    \end{tabular}
    \end{minipage}"""


def write_run_info(handle):
    """Write general information on the run."""
    run_info = """

Run information
---------------

Mutations at single nucleotide positions were called with a %4.1f%% frequency threshold and a minimum depth of %d reads.

""" % (100 * MIN_FRACTION, RAW_DEPTH_THRESHOLD)
    print(run_info, file=handle)


def write_md(org=None, mut_file='final.csv', subtype_file='subtype_evidence.csv'):
    """Write the markdown file."""
    logging.info('Writing report in markdown')
    rh = open('report.md', 'w')
    write_subtype_info(rh, subtype_file)

    if org == 'HIV':
        resistance_mutations = parse_drm()
        write_sierra_results(rh, mut_file)
    elif org == 'HCV':
        resistance_mutations = parse_ras()
        write_header_HCV(rh, resistance_mutations)
    elif org is None:
        md = """
No HIV/HCV read found
=====================
"""
        print(md, file=rh)
        return

    write_run_info(rh)
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

        # drms.drop(['commented', 'mut_y'], axis=1, inplace=True)
        drms.rename(columns={'mut_x': 'mut'}, inplace=True)

        cols = ['gene', 'pos', 'mut', 'freq', 'category']
        drms = drms[cols]
        drms.sort_values(by=['gene', 'pos', 'freq'], inplace=True,
                         ascending=[True, True, False])
        drms.to_csv('annotated_DRM.csv', index=False)

        print('Mutations detected', file=rh)
        print('==================\n', file=rh)
        print('Mutation lists from %s.\n' % mastercomments_version, file=rh)
        # write_header_HIV(rh, resistance_mutations)
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
            print('-' * len(gene) + '\n', file=rh)
            h1 = '| position | mutation | frequency [%] |      category      |'
            print(h1, file=rh)
            h2 = '|:{:-^8}:|:{:-^8}:|:{:-^13}:|:{:-^18}:|'
            print(h2.format('', '', '', ''), file=rh)
            for index, row in gene_muts.iterrows():
                # # same pos mut tuple must give same annota, probably redundant
                # assert group['category'].nunique() == 1, \
                #     group['category'].describe()
                int_freq = int(round(100 * row['freq'], 0))
                if int_freq == 0:
                    continue
                if 'major' in row['category'].lower() or 'nrti' in row['category'].lower():
                    mut_cat = '**%s**' % row['category']
                else:
                    mut_cat = row['category']
                print(
                    '|{: ^10}|{: ^10}|{: ^15}|{: ^20}|'.format(int(row['pos']),
                                                               row['mut'],
                                                               int_freq,
                                                               mut_cat),
                    file=rh)
            print('\n', file=rh)
    elif org == 'HCV':
        drms = pd.read_csv('merged_muts_drm_annotated.csv')
        drms = drms.fillna("unknown")
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
                if int_freq == 0:
                    continue
                category = row['CATEGORY']
                print(
                    '|{: ^10}|{: ^10}|{: ^15}|{: ^12}'.format(int(row['pos']), row['mut'], int_freq, category),
                    file=rh)
                del index
            print('\n', file=rh)
    write_ambig_score(rh)
    rh.close()


def convert_2_pdf(fastq=None, version='unknown'):
    """Convert markdown file to pdf with pandoc, filling sample and version info."""
    import shutil
    import re
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
    logging.info('Converting markdown to pdf with pandoc')
    pand_cml = 'pandoc --template=./template.tex report.md -o report.pdf'
    logging.debug(pand_cml)
    subprocess.call(pand_cml, shell=True)
    os.remove('./template.tex')
    os.remove('contact.tex')


def main(org=None, subtype_file=None, fastq='unknown', version='unknown'):
    """What the main does."""
    write_md(org=org, mut_file='final.csv', subtype_file=subtype_file)
    convert_2_pdf(fastq=fastq, version=version)


if __name__ == '__main__':
    #write_sierra_results(sys.stdout, 'final.csv')
    main(org=sys.argv[1], subtype_file='subtype_evidence.csv')
