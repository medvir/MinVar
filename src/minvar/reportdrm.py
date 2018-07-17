#!/usr/bin/env python3
"""

Make report with Drug Resistance Mutations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Parse the different output files and write a report in markdown, then convert it to pdf with pandoc.

This could be improved by writing pure functions that return strings of markdown formatted text and
letting the main orchestrate the composition of paragraphs.

"""
import sys
import os
import configparser
import csv
import logging
import shlex
import subprocess
from pkg_resources import resource_filename
from Bio import SeqIO, AlignIO

import pandas as pd

dn_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if __name__ == '__main__':
    if __package__ is None:
        os.sys.path.insert(1, dn_dir)
        mod = __import__('minvar')
        sys.modules["minvar"] = mod
        from common import MIN_FRACTION, RAW_DEPTH_THRESHOLD, drug_names, mastercomments_version
        from Alignment import needle_align
else:
    from .common import MIN_FRACTION, RAW_DEPTH_THRESHOLD, drug_names, mastercomments_version
    from .Alignment import needle_align

cons_B_file = resource_filename(__name__, 'db/HIV/consensus_B.fna')

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
    """Return drug resistance mutations listed in files db/HIV/masterComments*.txt."""
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
    """Return position of RAS listed in file db/HCV/all_mutations_position.csv."""
    ras_file = resource_filename(__name__, 'db/HCV/all_mutations_position.csv')
    ras_positions = pd.read_csv(ras_file)
    ras_positions['CATEGORY'] = 'RAS'
    return ras_positions


def write_subtype_info(handle, subtype_file=None):
    """If a subtype evidence file is present, write information on subtyping in markdown into ``handle``.

    Reads the support given to each subtype/genotype in ``subtype_file``. If the highest support
    is < 50%% then writes "Undetermined"

    :param subtype_file: path to csv file with subtype,reads columns
    :param handle: handle of markdown file being written
    """
    from operator import itemgetter

    try:
        save_freq = {}
        with open(subtype_file) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for mtype, freq in spamreader:
                int_freq = int(round(100 * float(freq), 0))
                if int_freq >= 1:
                    save_freq[mtype] = save_freq.get(mtype, 0) + int_freq
        top_gt, support = sorted(save_freq.items(), key=itemgetter(1), reverse=True)[0]
        if support >= 50:
            md_header = 'Subtype with the highest blast support: %s (%d %%)\n' % (top_gt, support)
            md_header += '-------------------------------------------------\n'
        else:
            md_header = 'Undetermined subtype: (blast support below 50%)\n'
            md_header += '---------------------------------------------\n'
    except FileNotFoundError:
        md_header = '\n'
    print(md_header, file=handle)


def write_sierra_results(handle, mut_file):
    """Run ``sierrapy patterns``, parse the results and writes them into ``handle``.

    Read mutations from ``mut_file`` and write only those above ``sierra_threshold`` them into a format that can be
    read with ``sierrapy patterns``.
    Then, read the resulting json and write it with the correct markdown formatting into ``handle``. This
    markdown also includes some latex to colour the table cells according to the resistance.

    This is only run on HIV samples. See `sierrapy page
    <https://github.com/hivdb/sierra-client/tree/master/python>`_ for more info.

    :param mut_file: path to csv file with gene,freq,wt,mut columns
    :param handle: handle of markdown file being written

    """
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
    print('Drug resistance interpretation', file=handle)
    print('==============================\n', file=handle)
    print('Stanford HIVdb version %s, pubdate %s.' % (version, pubdate), file=handle)
    print('Mutations below %d%% were not included.\n' % (100 * sierra_threshold), file=handle)

    for dr in pattern['drugResistance']:
        gene_name = dr['gene']['name']
        comments_set = set()
        print(gmap_inverse.get(gene_name, gene_name), file=handle)
        print('-' * len(gene_name) + '\n', file=handle)
        print('| class |         name         | score|      assessment      |              mutations               |',
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
    """Read ambiguous sequence, compute and write ambiguity score to handle.

    Read the sequence from ``cns_ambiguous.fasta``, align it to PRRT up to amminoacid 335, extract
    the aligned region and count the total and ambiguous nucleotides (*i.e.*, not ACGT).

    :param handle: handle of markdown file being written
    """
    print('\nAmbiguity score', file=handle)
    print('---------------\n', file=handle)
    amb_target = str(list(SeqIO.parse(cons_B_file, 'fasta'))[0].seq[168:1470])
    needle_align('cns_ambiguous.fasta', 'asis:%s'% amb_target, 'ambi_aln.fasta', go=10., ge=.5)
    alignment = AlignIO.read('ambi_aln.fasta', 'fasta')
    start = None
    i = 0
    while start is None:
        if '-' not in alignment[:, i]:
            start = i
        i += 1
    ambi_al, cons_al = str(alignment[0].seq), str(alignment[1].seq)
    stop = min(len(cons_al.rstrip('-')), len(ambi_al.rstrip('-')))
    target = ambi_al[start:stop].replace('-', '').replace('N', '')
    rlen = len(target)
    score = float(sum((1 for nt in target if nt not in set(['A', 'C', 'G', 'T']))))
    print('Score: %4.2f %% (%d of %d total nucleotides).\n\n' % (100 * score / rlen, score, rlen), file=handle)
    print('Region to compute ambiguity score is 1302 bp, reached high coverage on %d.\n\n' % rlen, file=handle)

# def write_header_HIV(handle, drms=None):
#     """Write header to a file in markdown format."""
#     md_header = """\
# Parsing mutations
# -----------------
#
# The list of annotated mutations was downloaded from HIVdb and includes:
#
# """
#     for gene in ['protease', 'RT', 'integrase']:
#         positions = set(drms[drms.gene == gene].pos.tolist())
#         md_header += '- %d positions on %s\n' % (len(positions), gene)
#     print(md_header, file=handle)


def write_header_HCV(handle, drms=None):
    """Write header with the count of known RAS to a file in markdown format (HCV only).

    :param handle: handle of markdown file being written
    :param drms: dataframe with known RAS mutation on each gene, colums: |gene|pos|
    """
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
    """OUTDATED: manual conversion of mutations extracted from HIVdb comment files.

    In older versions, comment files contained information with NOT XYZ to signify
    all aminoacids except X, Y, and Z. This function translated this into a list of
    valid ammino acids. It was too complicated to achieve this functionality with panda alone.
    """
    # read file with annotations, result of a merge detected with HIVdb annotations
    annotated = pd.read_csv(mer_file)
    # prepare output file
    genes = []
    positions = []
    wts = []
    muts = []
    freqs = []
    comments = []
    # iterate on annotations grouped by gene/pos
    for name, group in annotated.groupby(['gene', 'pos']):
        aas = group[pd.notna(group['mut_x'])]['mut_x'].tolist()
        if not aas:  # mutation not detected at this position
            continue
        mutation_found = False
        # mut_x is detected, mut_y is annotation, iterate on detected to see if it is among the annotated ones
        for idx, row in group.iterrows():
            detected_mutation = row['mut_x']
            if detected_mutation in str(row['mut_y']):
                mutation_found = True
                genes.append(name[0])
                positions.append(name[1])
                wts.append(row['wt'])
                freqs.append(row['freq'])
                muts.append(row['mut_x'])
                comments.append(row['comment'])
                break
        if not mutation_found:
            genes.append(name[0])
            positions.append(name[1])
            wts.append(row['wt'])
            freqs.append(row['freq'])
            muts.append(row['mut_x'])
            comments.append('unannoted')
    mdf = pd.DataFrame({'gene': genes, 'pos': positions, 'wt': wts, 'mut': muts, 'freq': freqs, 'comment': comments})
    print(len(mdf))
    return mdf


def write_contact_file(sample_id='unknown sample', version='unknown'):
    """Write tex file with contact information taken from ini file and sample_id in footer.

    :param sample_id: this will be written in the footer, left
    :param version: MinVar version to be written in the footer, right
    """
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
    oh.write(r'%s\\' % contact_dict['unit'] + '\n')
    oh.write(r'phone: %s\\' % contact_dict['phone'] + '\n')
    oh.write(r'fax: %s\\' % contact_dict['fax'] + '\n')
    oh.write(r'email: \href{mailto:%s}{%s}\\' % (contact_dict['email'], contact_dict['email']) + '\n')
    oh.write(r'\end{tabular}' + '\n')
    oh.write(r'\end{minipage}' + '\n')
    oh.write(r'\vspace{1cm}' + '\n')

    # r"""\fancyfoot[L]{John Smith \quad - \quad \href{mailto:john@smith.com}{{\small john@smith.com}}}
    #
    # \begin{minipage}{0.5\textwidth}
    #   \includegraphics[width=2in]{uzh_logo_e_pos}% Logo
    # \end{minipage}
    # \hfill
    # \begin{minipage}{0.45\textwidth}
    # \begin{tabular}{@{}r@{}}
    #   %\DTMnow \\[\normalbaselineskip]
    #   Institute of Medical Virology \\
    #   Telephone + 41 44 63 42653  \\
    #   Fax + 41 44 63 44967 \\
    #   Email \href{mailto:med-virology@virology.uzh.ch}{med-virology@virology.uzh.ch}
    # \end{tabular}
    # \end{minipage}"""


def write_run_info(handle):
    """Write general information on the run. Possible future expansion."""
    run_info = """

Run information
---------------

Mutations at single nucleotide positions were called with a %4.1f%% frequency
threshold and a minimum depth of %d reads.

""" % (100 * MIN_FRACTION, RAW_DEPTH_THRESHOLD)
    print(run_info, file=handle)


def write_md(org=None, mut_file='final.csv', subtype_file='subtype_evidence.csv', sample_id=''):
    """Orchestrate the writing of the markdown file.

    :param org: HIV or HCV
    :param mut_file: csv file where ammino acid mutations are written with columns gene,wt,pos,mut,freq
    :param subtype_file: csv file with subtype,support columns
    :param sample_id: will be written in the title of the report
    """
    logging.info('Writing report in markdown')
    rh = open('report.md', 'w')
    print('Sample %s: drug resistance mutations report' % sample_id, file=rh)
    print('============================================\n', file=rh)
    # print('Drug resistance mutations detected by NGS sequencing', file=rh)
    # print('====================================================\n', file=rh)

    write_subtype_info(rh, subtype_file)

    if org == 'HIV':
        resistance_mutations = parse_drm()
        write_sierra_results(rh, mut_file)
    elif org == 'HCV':
        resistance_mutations = parse_ras()
        # write_header_HCV(rh, resistance_mutations)
    elif org is None:
        md = """
No HIV/HCV read found
=====================
"""
        print(md, file=rh)
        return

    logging.info('Reading mutations from %s', mut_file)
    mutation_detected = pd.read_csv(mut_file)
    logging.info('shape: %s', mutation_detected.shape)
    logging.info('Parsed DRM from database, shape: %s', str(resistance_mutations.shape))
    mpd = pd.merge(mutation_detected, resistance_mutations, how='left',
                   on=['gene', 'pos'])
    mpd.to_csv(path_or_buf='merged_muts_drm_annotated.csv', index=False)
    logging.info('Shape of raw merged is: %s', str(mpd.shape))

    print('Mutations detected', file=rh)
    print('==================\n', file=rh)
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

        print('Mutation lists from %s.\n' % mastercomments_version, file=rh)
        write_run_info(rh)
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
                    '|{: ^10}|{: ^10}|{: ^15}|{: ^20}|'.format(int(row['pos']), row['mut'], int_freq, mut_cat),
                    file=rh)
            print('\n', file=rh)
        write_ambig_score(rh)
    elif org == 'HCV':
        write_header_HCV(rh, resistance_mutations)
        write_run_info(rh)
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
    rh.close()


def convert_2_pdf(sample_id='', version='unknown'):
    """Convert markdown file to pdf with pandoc, filling sample and version info.

    :param sample_id: sample name
    :param version: MinVar version
    """
    import shutil

    # copy template to current directory
    tmpl_file = resource_filename(__name__, 'db/template.tex')
    shutil.copy(tmpl_file, os.getcwd())
    # write contact.tex to be included in template
    write_contact_file(sample_id=sample_id, version=version)
    logging.info('Converting markdown to pdf with pandoc')
    pand_cml = 'pandoc --template=./template.tex report.md -o report.pdf'
    logging.debug(pand_cml)
    subprocess.call(pand_cml, shell=True)
    os.remove('./template.tex')
    os.remove('contact.tex')


def main(org=None, subtype_file=None, fastq='unknown', version='unknown'):
    """What the main does."""
    import re

    fastq_base = os.path.basename(fastq)
    try:
        sample_id = re.search(r'(.*)_S\d*', fastq_base).group(1)
    except AttributeError:
        sample_id = 'unknown'
    write_md(org=org, mut_file='final.csv', subtype_file=subtype_file, sample_id=sample_id)
    convert_2_pdf(sample_id=sample_id, version=version)


if __name__ == '__main__':
    # write_sierra_results(sys.stdout, 'final.csv')
    x = parse_merged('merged_muts_drm_annotated.csv')
    x.to_csv('x.csv')
    #main(org=sys.argv[1], subtype_file='subtype_evidence.csv')
