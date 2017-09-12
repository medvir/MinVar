#!/usr/bin/env python3

import os
import sys
import unittest
import tempfile

minvar_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(1, minvar_dir)
mod = __import__('minvar')
sys.modules["minvar"] = mod

from minvar.prepare import compute_min_len, filter_reads, iterate_consensus, find_subtype


class Test_compute_min_len(unittest.TestCase):

    def setUp(self):
        self.hiv_fastq = '/Users/ozagordi/hiv_small.fastq.gz'

    def test_minlen_one(self):
        m = compute_min_len(self.hiv_fastq)
        self.assertGreater(m, 49)


class Test_iterate_consensus(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.gettempdir()
        self.hiv_fastq = '/Users/ozagordi/hiv_small.fastq.gz'
        self.hiv_ref = '/Users/ozagordi/References/HIV-HXB2-pol.fasta'
        #self.hcv_ref = '/Users/ozagordi/References/HIV-HXB2-pol.fasta'

    # def test_filter_only(self):
    #     os.chdir(self.tmpdir)
    #     hq = filter_reads(self.hiv_fastq, 1000, 50)
    #     self.assertIsInstance(hq, str)

    def test_iterate_consensus(self):
        os.chdir(self.tmpdir)
        cons = iterate_consensus(self.hiv_fastq, self.hiv_ref)
        self.assertIsInstance(cons, str)


class Test_find_subtype(unittest.TestCase):

    def setUp(self):
        self.tmpdir = '/tmp/ppp' #tempfile.gettempdir()
        self.hiv_fastq = '/Users/ozagordi/hiv_small.fastq.gz'
        self.hcv_fastq = '/Users/ozagordi/hcv_small.fastq.gz'

    def test_hiv(self):
        os.chdir(self.tmpdir)
        org, bs, sf = find_subtype(self.hiv_fastq)
        self.assertEqual(org, 'HIV')
        self.assertEqual(bs, 'CONSENSUS_B')

    def test_hcv(self):
        os.chdir(self.tmpdir)
        org, bs, sf = find_subtype(self.hcv_fastq)
        self.assertEqual(org, 'HCV')
        self.assertEqual(bs, 'M62321.1')

if __name__ == '__main__':
    unittest.main()
