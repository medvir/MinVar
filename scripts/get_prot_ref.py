#!/usr/bin/env python3
'''Sequences in refs are full genome from different genotypes.
This looks for ref_id in refs, translates and aligns to a protein from
h77 in order to extract the one for the requested genotype
'''
import os
import sys
import shlex
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# manipulate path to import functions
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(1, parent_dir)
mod = __import__('minvar')
sys.modules["minvar"] = mod
from minvar.common import hcv_map

def extract_protein(ref_id, prot):
    '''Workhorse for the script, takes a single accession number and protein
    and returns the corresponding translated region by aligning against the
    protein from H77 isolate
    '''

    try:
        os.remove('aln.fasta')
    except FileNotFoundError:
        pass
    try:
        os.remove('tmp.fasta')
    except FileNotFoundError:
        pass
    ref_file = os.path.join(parent_dir, 'minvar/db/HCV/HCV_ref.fasta')
    refs = list(SeqIO.parse(ref_file, 'fasta'))

    raw_file = os.path.join(parent_dir, 'minvar/db/HCV/%s_h77.faa' % prot)
    #raw_file = list(SeqIO.parse(raw_file, 'fasta'))[0].seq

    for r in refs:
        if r.id.startswith(ref_id):
            ref = r
        else:
            continue

    frame = []
    stops = 10000
    for f in range(3):
        frame.append(ref.seq[f:].translate())
        stops_here = str(frame[f]).count('*')
        if stops_here < stops:
            stops = stops_here
            best_frame = f

    orf = str(frame[best_frame])
    sr = SeqRecord(Seq(orf), id='best_orf', description='')
    SeqIO.write(sr, 'tmp.fasta', 'fasta')
    cml = 'needle -asequence %s -bsequence tmp.fasta' % shlex.quote(raw_file)
    cml += ' -aformat3 fasta -auto -outfile aln.fasta'

    subprocess.call(shlex.split(cml))
    aln = AlignIO.read('aln.fasta', 'fasta')
    base = str(aln[0, :].seq)
    target = str(aln[1, :].seq)
    start = len(base) - len(base.lstrip('-'))
    stop = len(base.rstrip('-'))
    dist = sum((base[i] != target[i] for i in range(start, stop)))
    dist = float(dist) / (stop - start)
    if dist > 0.0:
        print(round(dist, 3))
        print(ref_id, prot)
        print(target[start:stop])
        print(base[start:stop])
    return target[start:stop]


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

h77_file = os.path.join(parent_dir, 'minvar/db/HCV/H77_polyprotein.faa')
h77_poly = list(SeqIO.parse(h77_file, 'fasta'))[0]
for k, v in h77_locations.items():
    sta, sto = v
    s_id = '%s_h77' % k
    sr1 = SeqRecord(h77_poly[sta:sto].seq, id=s_id, description='')
    SeqIO.write(sr1, s_id + '.faa', 'fasta')

s_ids = ['M62321', 'M58335', 'D00944', 'D10988', 'D17763', 'Y11604']
target_prots = ['ns3', 'ns4a', 'ns4b', 'ns5a', 'ns5b']

for target_prot in target_prots:
    for s_id in s_ids:
        xt_prot = extract_protein(s_id, target_prot)
        gt = hcv_map[s_id]
        sr2 = SeqRecord(Seq(xt_prot), id='genotype_%s_%s' % (gt, target_prot),
                        description='')
        SeqIO.write(sr2, '%s_%s.faa' % (target_prot, gt), 'fasta')
