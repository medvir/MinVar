#!/usr/bin/env python3

import os
import tempfile

import pytest
from Bio import SeqIO

from src.minvar.prepare import (compute_min_len, find_subtype,
                                iterate_consensus, pad_consensus)

minvar_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cbfile = os.path.join(minvar_dir, 'src/minvar/db/HIV/consensus_B.fna')
hiv_middle = list(SeqIO.parse(cbfile, 'fasta'))[0]
rec_hcv_file = os.path.join(minvar_dir,
                            'src/minvar/db/HCV/recomb_references.fasta')
hcv_middle = list(SeqIO.parse(rec_hcv_file, 'fasta'))[0]


# @pytest.mark.skip(reason="no way of currently testing this")
def test_HIV_length():
    # a sequence from the middle of the reference
    m = pad_consensus(hiv_middle[100:500], 'HIV', 'whatever')
    assert len(m) == 3012
    # from the beginning
    m = pad_consensus(hiv_middle[:500], 'HIV', 'whatever')
    assert len(m) == 3012
    # from the end
    m = pad_consensus(hiv_middle[-500:], 'HIV', 'whatever')
    assert len(m) == 3012
    # longer on the left
    lotl = 'AGACTAGCCGATCAGCATCAGCA' + hiv_middle[:100]
    m = pad_consensus(lotl, 'HIV', 'whatever')
    assert len(m) == 3012
    # longer on the right
    lotl = hiv_middle[-100:] + 'AAGCGCATCGACATCAGCA'
    m = pad_consensus(lotl, 'HIV', 'whatever')
    assert len(m) == 3012


@pytest.mark.skip(reason="no way of currently testing this")
def test_HCV_recomb_length():
    m = pad_consensus(hcv_middle, 'HCV', 'RF1_2k/1b')
    assert len(m) == 9357

    m = pad_consensus(hcv_middle[:256], 'HCV', 'RF1_2k/1b')
    assert len(m) == 9357

    m = pad_consensus(hcv_middle[-256:], 'HCV', 'RF1_2k/1b')
    assert len(m) == 9357


def test_compute_min_len():
    """Test minimum length"""
    hiv_fastq = '/Users/ozagordi/hiv_small.fastq.gz'
    m = compute_min_len(hiv_fastq)
    assert m > 49


@pytest.mark.skip(reason="skipping iterate_consensus")
def test_iterate_consensus():

    tmpdir = tempfile.gettempdir()
    hiv_fastq = '/Users/ozagordi/hiv_small.fastq.gz'
    hiv_ref = '/Users/ozagordi/References/HIV-HXB2-pol.fasta'
    # hcv_ref = '/Users/ozagordi/References/HIV-HXB2-pol.fasta'

    # def test_filter_only(self):
    #     os.chdir(tmpdir)
    #     hq = filter_reads(hiv_fastq, 1000, 50)
    #     assertIsInstance(hq, str)

    os.chdir(tmpdir)
    cons = iterate_consensus(hiv_fastq, hiv_ref)
    assert isinstance(cons, str)


@pytest.mark.skip("skipping find_subtype")
def test_find_subtype():

    tmpdir = tempfile.gettempdir()
    hiv_fastq = '/Users/ozagordi/hiv_small.fastq.gz'
    hcv_fastq = '/Users/ozagordi/hcv_small.fastq.gz'

    # test hiv
    os.chdir(tmpdir)
    org, bs, sf = find_subtype(hiv_fastq)
    assert org == 'HIV'
    assert bs == 'CONSENSUS_B'

    # test hcv
    os.chdir(tmpdir)
    org, bs, sf = find_subtype(hcv_fastq)
    assert org == 'HCV'
    assert bs == 'M62321.1'
