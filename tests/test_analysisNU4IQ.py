import argparse
import tempfile
import unittest

import numpy

from analysisNU4IQ import main


class TestRegressionCsvOutput(unittest.TestCase):

	def test_regressionCsvOutput(self):
		expectedData = numpy.loadtxt('data/results.csv', delimiter=',', skiprows=2, usecols=range(1, 19))
		with tempfile.NamedTemporaryFile() as tmpFile:
			args = argparse.Namespace(iFile='data/rawValues.tsv', oFile=tmpFile.name, imName=None)
			main(args)
			obtainedData = numpy.loadtxt(tmpFile.name, delimiter=',', skiprows=2, usecols=range(1, 19))
			self.assertTrue(numpy.allclose(expectedData, obtainedData))


if __name__ == '__main__':
	unittest.main()
