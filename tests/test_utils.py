#
# tests/test_utils.py
#   - copy and/or modify together with dsdobjects/utils.py
#
import logging
logger = logging.getLogger('dsdobjects')
logger.setLevel(logging.INFO)
import unittest

from dsdobjects.utils import (flint, convert_units)

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class TestUtils(unittest.TestCase):
    def test_convert_units(self):
        assert convert_units(3, 'M', 'nM') == 3_000_000_000
        assert convert_units(2.4, 'M', 'nM') == 2_400_000_000
        assert convert_units(2.4, 'nM', 'mM') == 2.4e-6
        assert convert_units(10, 'nM', 'mM') == 1e-5

if __name__ == '__main__':
    unittest.main()
