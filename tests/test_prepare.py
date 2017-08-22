#!/usr/bin/env python3

import os
import sys
import unittest
import tempfile

minvar_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(1, minvar_dir)
mod = __import__('minvar')
sys.modules["minvar"] = mod

from minvar.prepare import compute_min_len, filter_reads, iterate_consensus


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
        cons = iterate_consensus(self.fastq_sample)
        self.assertIsInstance(cons, str)


if __name__ == '__main__':
    print('xyz')
    unittest.main()
