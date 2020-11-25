#
# tests/test_base_classes.py
#   - copy and/or modify together with dsdobjects/base_classes.py
#
import logging
logger = logging.getLogger('dsdobjects')
logger.setLevel(logging.INFO)
import unittest

import gc
from dsdobjects import (SingletonError,
                        show_singletons,
                        clear_singletons,
                        DomainS)
from dsdobjects.objectio import read_pil_line, set_io_objects, clear_io_objects

SKIP = True # Below test imports matplotlib, that will ruin later tests.

class MyDomain(DomainS):
    pass

def initdomain1():
    return read_pil_line('length a = 15')
    
def initdomain2():
    return read_pil_line('length a = 10')

class TestMemoryLeaks(unittest.TestCase):
    def setUp(self):
        set_io_objects(D = MyDomain)

    def tearDown(self):
        clear_io_objects()

    def test_memory_ok(self):
        for s in show_singletons(MyDomain):
            print(s)
        for _ in range(2):
            a = initdomain1()
        del a
        for _ in range(2):
            a = initdomain2()
        del a
        gc.collect()
        for s in show_singletons(MyDomain):
            print(s)

@unittest.skipIf(SKIP, "skipping tests.")
class TestMemoryLeaks2(unittest.TestCase):
    def setUp(self):
        set_io_objects(D = MyDomain)
        self.SKIP = False
        try:
            import matplotlib.pyplot
        except ImportError:
            self.SKIP = True

    def tearDown(self):
        clear_io_objects()
        clear_singletons(MyDomain)

    def test_memory_leak(self):
        if self.SKIP:
            return
        with self.assertRaises(SingletonError):
            for _ in range(2):
                a = initdomain1()
            del a
            for _ in range(2):
                a = initdomain2()
            del a

    def test_memory_leak_fix(self):
        for _ in range(2):
            clear_singletons(MyDomain)
            a = initdomain1()
        del a

        for _ in range(5):
            clear_singletons(MyDomain)
            a = initdomain2()
        del a

if __name__ == '__main__':
    unittest.main()

