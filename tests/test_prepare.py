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
        self.tmpdir = tempfile.gettempdir()
        self.fastq_sample = '/Users/ozagordi/hiv_small.fastq.gz'

    def test_minlen_one(self):
        m = compute_min_len(self.fastq_sample)
        self.assertGreater(m, 49)

    def test_filter_only(self):
        os.chdir(self.tmpdir)
        hq = filter_reads(self.fastq_sample, 1000, 50)
        self.assertIsInstance(hq, str)

    def test_iterate_consensus(self):
        os.chdir(self.tmpdir)
        cons = iterate_consensus(self.fastq_sample)
        self.assertIsInstance(cons, str)


class Test_find_subtype(unittest.TestCase):

    def setUp(self):
        self.tmpdir = '/tmp/ppp' #tempfile.gettempdir()
        self.hiv_fastq = '/Users/ozagordi/hiv_small.fastq.gz'
        self.hcv_fastq = '/Users/ozagordi/hcv_small.fastq.gz'

    def test_hiv(self):
        os.chdir(self.tmpdir)
        org, sf = find_subtype(self.hiv_fastq)
        self.assertEqual(org, 'HIV')

    def test_hcv(self):
        os.chdir(self.tmpdir)
        org, sf = find_subtype(self.hcv_fastq)
        self.assertEqual(org, 'HCV')


if __name__ == '__main__':
    unittest.main()
