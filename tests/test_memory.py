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
import dsdobjects.objectio as oio

class MyDomain(DomainS):
    pass

oio.Domain = MyDomain

def initdomain1():
    return oio.read_pil_line('length a = 15')
    
def initdomain2():
    return oio.read_pil_line('length a = 10')

class TestMemoryLeaks(unittest.TestCase):
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

class TestMemoryLeaks2(unittest.TestCase):
    def test_memory_leak(self):
        with self.assertRaises(SingletonError):
            import matplotlib.pyplot
            for _ in range(2):
                a = initdomain1()
            del a
            for _ in range(2):
                a = initdomain2()
            del a

    def test_memory_leak_fix(self):
        import matplotlib.pyplot
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

