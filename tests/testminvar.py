import unittest

class TestMakeCons(unittest.TestCase):

	def setUp(self):
		pass

	def testfailintentionally(self):
		self.assertEqual(1 == 2)


if __name__ == '__main__':
	unittest.main()
