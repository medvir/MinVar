import unittest

class TestMakeCons(unittest.TestCase):

	def setUp(self):
		pass

	def testfailintentionally(self):
		self.assertEqual(1 == 2)

	def testbasicpass(self):
		self.assertEqual('abc' == 'abc')


if __name__ == '__main__':
	unittest.main()
