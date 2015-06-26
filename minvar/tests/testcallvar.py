#!/usr/local/bin/python3.4
import os
import sys
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

import pandas as pd

minvar_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.sys.path.insert(1, minvar_dir)
mod = __import__('minvar')
sys.modules["minvar"] = mod

from ReportDRM.src import Alignment
from minvar import callvar

class TestCallVar(unittest.TestCase):

	def setUp(self):
		# two silent mutations w.r.t. consensus B
		self.ref_file_1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mock_ref_1.fasta')
		self.ref_1 = SeqRecord(Seq('TTATCTATCAGTACATGGATGATTTGTATGTAGGATCTGACTTAGAAATTGGGCAGCATAGAACAAAAATAGAGGAACTGAGACAACATCTGTTGAGGTGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTT', generic_dna), id='mock_1')
		SeqIO.write(self.ref_1, self.ref_file_1, 'fasta')
		# two mutation w.r.t. consensus B, one silent
		self.ref_file_2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mock_ref_2.fasta')
		self.ref_2 = SeqRecord(Seq('TTATCTATCAGTACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGACAACATCTGTTGAGGTGCGGATTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTT', generic_dna), id='mock_2')
		SeqIO.write(self.ref_2, self.ref_file_2, 'fasta')
		# frame is 1 here
		self.ref_file_3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mock_ref_3.fasta')
		self.ref_3 = SeqRecord(Seq('TTTTTTAGGGAAGATCTGGCCTTCCCACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAAGAGACAACAACTCCCTCTCAGAAGCAGGAGCCGATAGACAAGGA', generic_dna), id='mock_3')
		SeqIO.write(self.ref_3, self.ref_file_3, 'fasta')

		self.test_vcf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vcf_test.vcf')
		assert os.path.exists(self.test_vcf)

	def test_parse_cons_mutations(self):
		cm = callvar.parse_cons_mutations(self.ref_1)
		cm = cm[cm.freq > 0.0]
		self.assertEqual(len(cm), 2, cm)
		#self.fail("Test to be written")

	def test_parsevar(self):
		out = callvar.parsevar(self.test_vcf, self.ref_file_1)
		self.assertEqual(len(out), 2, out)
		out_0 = out.iloc[0]
		out_1 = out.iloc[1]
		self.assertEqual(out_0['wt'], 'T', out)
		self.assertEqual(out_0['pos'], 4, out)
		self.assertEqual(out_0['mut'], 'A', out)
		self.assertEqual(out_1['wt'], 'G', out)
		self.assertEqual(out_1['pos'], 11, out)
		self.assertEqual(out_1['mut'], 'C', out)

	def test_find_frame(self):
		frame, nt_framed, aa_framed = callvar.find_frame(self.ref_1)
		self.assertEqual(frame, 3)
		self.assertEqual(len(nt_framed) % 3, 0)
		self.assertEqual(nt_framed[-1], 'T', nt_framed)
		self.assertEqual(aa_framed, 'IYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFL')

	def test_merge_mutations(self):
		cm = callvar.parse_cons_mutations(self.ref_1)
		'''
		A, 11, G
		A, 50, T
		'''
		vm = pd.DataFrame([['T', 1, 'A', 0.3],  # include new mutation
						   ['G', 11, 'T', 0.2],  # already mutated in sample cons
			               ['T', 50, 'A', 0.1]], # mutates back to cons B
			              columns=['wt', 'pos', 'mut', 'freq'])
		m_out = callvar.merge_mutations(cm, vm)
		# two mutations in consensus, three in vcf but the one in pos 101 is
		# reverting to consensus B
		self.assertEqual(len(m_out), 4, m_out)
		freqs = [0.3, 0.8, 0.2, 0.9]
		muts = ['A', 'G', 'T', 'T']
		for i in range(4):
			out_here = m_out.iloc[i]
			with self.subTest(mut=i):
				self.assertEqual(out_here.freq, freqs[i])
				self.assertEqual(out_here.mut, muts[i])

	def test_annotate_mutations(self):
		import pandas as pd
		muts = pd.DataFrame(columns=['gene', 'pos', 'wt', 'mut', 'freq'])

		# no mutation first
		out = callvar.annotate_mutations(muts, self.ref_3)
		self.assertEqual(len(out), 0)

		# add one mutation
		muts = muts.append({'wt': 'T', 'pos': 6, 'mut': 'A', 'freq': 0.2}, ignore_index=True)
		out = callvar.annotate_mutations(muts, self.ref_3)
		self.assertEqual(len(out), 1, out)
		out_0 = out.iloc[0]
		self.assertEqual(out_0.pos, 2)
		self.assertEqual(out_0.gene, 'GagPolTF')
		self.assertEqual(out_0.wt, 'F')
		self.assertEqual(out_0.mut, 'L')
		self.assertEqual(out_0.freq, 0.2)

	def tearDown(self):
		os.remove(self.ref_file_1)
		os.remove(self.ref_file_2)
