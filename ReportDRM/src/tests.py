import reportdrm
import random
import unittest

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

class TestReportDRMFunctions(unittest.TestCase):

    def setUp(self):
        self.haplos = []
        # modify first residue in protease P -> Q
        seq = Seq('QQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD\
QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF', generic_dna)
        self.haplos.append(SeqRecord(seq, id='hap1', description='hap1freq=0.123'))
        # delet the second residue, change the third to F,
        seq = Seq('QFTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD\
QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF', generic_dna)
        self.haplos.append(SeqRecord(seq, id='hap2', description='hap2 freq=0.234'))

    def test_parse_mutations(self):
        pms = reportdrm.parse_mutations(self.haplos)

        self.assertTrue(pms.shape[0] >= len(self.haplos))
    def test_write_more_tests(self):
        self.assertEqual('Write more tests', '===')

if __name__ == '__main__':
    unittest.main()