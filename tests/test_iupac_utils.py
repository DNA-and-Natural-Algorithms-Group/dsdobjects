#
# tests/test_iupac_utils.py
#   - copy and/or modify together with dsdobjects/iupac_utils.py
#
import logging
logger = logging.getLogger('dsdobjects')
logger.setLevel(logging.INFO)
import unittest

from dsdobjects.iupac_utils import complement

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class TestConstraints(unittest.TestCase):
    def test_complement(self):
        sequence = 'AAAA'
        assert 'TTTT' == complement(sequence)

if __name__ == '__main__':
    unittest.main()
