#!/usr/bin/env python3
"""Script used to automatically obtain several references.

Sequences in refs are full genome from different genotypes.
This looks for ref_id in refs, translates and aligns to a protein from
h77 in order to extract the one for the requested genotype
"""
import os
import os.path
import sys
import shlex
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# manipulate path to import functions
if __name__ == '__main__':
    if __package__ is None:
        os.sys.path.insert(1, parent_dir)
        mod = __import__('minvar')
        sys.modules["minvar"] = mod
from common import hcv_map, acc_numbers

# zero-based coordinates on AAB67038.1 polyprotein from H77 isolate
h77_locations = {
    'c': (0, 191),
    'e1': (191, 383),
    'e2': (383, 746),
    'p2': (746, 809),
    'ns2': (809, 1026),
    'ns3': (1026, 1657),
    'ns4a': (1657, 1711),
    'ns4b': (1711, 1972),
    'ns5a': (1972, 2420),
    'ns5b': (2420, 3011)
}

target_prots = ['ns3', 'ns4a', 'ns4b', 'ns5a', 'ns5b']


def extract_protein(polyprotein_file, prot):
    """Obtain specific proteins by aligning against H77.

    Workhorse for the script, takes a single accession number and protein
    and returns the corresponding translated region by aligning against the
    protein from H77 isolate
    """
    cml2 = 'needle -asequence %s_h77.faa -bsequence %s' % (prot,
                                                           polyprotein_file)
    cml2 += ' -aformat3 fasta -auto -outfile aln.fasta'
    subprocess.call(shlex.split(cml2))
    aln = AlignIO.read('aln.fasta', 'fasta')
    os.remove('aln.fasta')
    base = str(aln[0, :].seq)
    target = str(aln[1, :].seq)
    start = len(base) - len(base.lstrip('-'))
    stop = len(base.rstrip('-'))
    dist = sum((base[i] != target[i] for i in range(start, stop)))
    dist = float(dist) / (stop - start)
    if dist > 0.0:
        print(round(dist, 3))
        print(target[start:stop])
        print(base[start:stop])
    return target[start:stop]


def split_h77():
    """Make all H77 proteins."""
    print('Splitting H77')
    h77_aa_file = os.path.join(parent_dir, 'db/HCV/H77_polyprotein.faa')
    if not os.path.exists(h77_aa_file):
        cml = 'efetch -db nuccore -format fasta -id AAB66324.1'
        with open(h77_aa_file, 'w') as f:
            subprocess.call(shlex.split(cml), stdout=f)
    h77_poly = list(SeqIO.parse(h77_aa_file, 'fasta'))[0]
    for k, v in h77_locations.items():
        sta, sto = v
        s_id = '%s_h77.faa' % k
        fpath = os.path.join(parent_dir, 'db/HCV/%s' % s_id)
        sr1 = SeqRecord(h77_poly[sta:sto].seq, id=s_id, description='')
        SeqIO.write(sr1, fpath, 'fasta')


def make_subtype_proteins(poly_files):
    """Extract proteins from references aligning to H77 proteins."""
    #s_ids = [d[0] for d in acc_numbers.items()]
    for target_prot in target_prots:
        for poly_file in poly_files:
            xt_prot = extract_protein(poly_file, target_prot)
            print(poly_file)
            gt, ac = os.path.split(poly_file)[1].split('_')[:2]
            sr2 = SeqRecord(Seq(xt_prot),
                            id='genotype_%s_%s_%s' % (gt, target_prot, ac),
                            description='')
            fpath = os.path.join(parent_dir,
                                 'db/HCV/%s_%s.faa' % (target_prot, gt))
            SeqIO.write(sr2, fpath, 'fasta')


def get_all_polyproteins():
    """Extract polyprotein from all references."""
    poly_files = []
    for s, v in hcv_map.items():
        cml = 'efetch -db nuccore -format gp -id %s -mode xml' % s
        cml += ' | xtract -pattern GBSeq -element GBSeq_locus -division GBFeature'
        cml += ' -block GBFeature_quals -section GBQualifier -if GBQualifier_name'
        cml += ' -equals translation -element GBQualifier_value'
        proc = subprocess.Popen(cml, shell=True, stdout=subprocess.PIPE,
                                universal_newlines=True)
        while True:
            line = proc.stdout.readline()
            if line:
                acn, sequence = line.split()
                if len(sequence) < 1000:
                    break
                print('Writing %s %s, length: %d' % (v, acn, len(sequence)))
                sr = SeqRecord(Seq(sequence), id=acn, description='polyprotein')
                fn = '%s_%s.faa' % (v.replace('/', '-'), acn.split('.')[0])
                SeqIO.write(sr, fn, 'fasta')
                poly_files.append(fn)
            else:
                break
    return poly_files


def get_references():
    """Download HCV references defined in common.py plus H77."""
    from pprint import pprint
    genotype = {}
    sub_list = []
    rec_list = []
    for k, v in acc_numbers.items():
        for n in v:
            genotype[n] = k
            if len(k) == 2:  # non recombinant subtypes: 1a, 1b, 1d
                sub_list.append(n)
            elif len(k) > 6:
                rec_list.append(n)
            else:
                print("Don't know what to do!")
    pprint(genotype)

    sub_path = os.path.join(parent_dir, 'db/HCV/subtype_references.fasta')
    if not os.path.exists(sub_path):
        cml = 'efetch -db nuccore -format fasta -id \"%s\"' % ','.join(sub_list)
        with open(sub_path, 'w') as f:
            subprocess.call(shlex.split(cml), stdout=f)

    rec_path = os.path.join(parent_dir, 'db/HCV/recomb_references.fasta')
    if not os.path.exists(rec_path):
        cml = 'efetch -db nuccore -format fasta -id \"%s\"' % ','.join(rec_list)
        with open(rec_path, 'w') as f:
            subprocess.call(shlex.split(cml), stdout=f)

    h77_nt_file = os.path.join(parent_dir, 'db/HCV/H77_cds.fna')
    if not os.path.exists(h77_nt_file):
        cml = 'efetch -db nuccore -format fasta -id AF009606.1'
        cml += ' -seq_start 342 -seq_stop 9374'
        with open(h77_nt_file, 'w') as f:
            subprocess.call(shlex.split(cml), stdout=f)

if __name__ == '__main__':
    split_h77()
    get_references()
    all_poly = get_all_polyproteins()
    make_subtype_proteins(all_poly)
