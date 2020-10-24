#
# tests/test_auto_classes.py
#   - copy and/or modify together with dsdobjects/auto_classes.py
#
import logging
logger = logging.getLogger('dsdobjects')
logger.setLevel(logging.INFO)
import unittest

from dsdobjects import ObjectInitError, SingletonError, clear_singletons
from dsdobjects.auto_classes import DomainA, ComplexA

SKIP = False
SKIP_TOXIC = True

@unittest.skipIf(SKIP, "skipping tests.")
class TestAutomaticDomain(unittest.TestCase):
    def tearDown(self):
        clear_singletons(DomainA)

    @unittest.skipIf(SKIP_TOXIC, "skipping toxic tests.")
    def test_garbage_collection(self):
        # This test must fail to see if garbage collection works.
        # Spoiler: it doesn't that is why we need tearDown()
        a = DomainA('a', 10)
        assert False

    def test_initialization_01(self):
        # Not enough parameters to initialize
        with self.assertRaises(SingletonError):
            a = DomainA('a')
        a = DomainA('a', 10)
        assert a is DomainA('a')
        assert DomainA('a', 10) is ~DomainA('a*', 10)

    def test_initialization_02(self):
        DomainA.ID = 1
        a = DomainA(dtype = 'short')
        b = DomainA(prefix = 'a', length = 10)
        c = DomainA(prefix = 'a', length = 10)
        assert (a.name, a.length) == ('d1', 5)
        assert (b.name, b.length) == ('a2', 10)
        assert (c.name, c.length) == ('a3', 10)

    def test_immutable_forms(self):
        a = DomainA('a', 15)
        b = DomainA('b', 15)

        # Forbidden change of immutable attributes.
        with self.assertRaises(SingletonError):
            a.length = 10
        assert len(a) == 15

        with self.assertRaises(SingletonError):
            a.name = 'b'
        assert a.name == 'a'

        with self.assertRaises(SingletonError):
            a = DomainA('b', 10)

@unittest.skipIf(SKIP, "skipping tests.")
class TestAutomaticComplex(unittest.TestCase):
    def setUp(self):
        self.d1 = DomainA(length = 20)
        self.d2 = DomainA(length = 5)

    def tearDown(self):
        clear_singletons(DomainA)
        clear_singletons(ComplexA)

    @unittest.skipIf(SKIP_TOXIC, "skipping toxic tests.")
    def test_garbage_collection(self):
        # This test must fail to see if garbage collection works.
        # Spoiler: it doesn't that is why we need tearDown()
        a = ComplexA(list('ABCDEFG'), list('.......')) 
        assert False

    def test_initialization_01(self):
        # Not enough parameters to initialize
        ComplexA.ID = 1
        with self.assertRaises(SingletonError):
            a = ComplexA(None, None, name = 'hello')
        a = ComplexA(list('ABCDEFG'), list('.......')) 
        assert a.kernel_string == 'A B C D E F G'
        assert str(a) == 'cplx1'
        b = ComplexA(list('A'), list('.'))
        assert b.kernel_string == 'A'
        assert str(b) == 'cplx2'
        assert a is ComplexA(None, None, name = 'cplx1')
    
    def test_initialization_02(self):
        d1, d2 = self.d1, self.d2
        a = ComplexA([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))
        assert a is ComplexA([d2, ~d1, '+', d1, d1, ~d2], list('((+).)'), name = a.name)
        assert a is ComplexA(None, None, name = a.name)

    def test_exceptions(self):
        d1, d2 = self.d1, self.d2
        a = ComplexA([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))
        try:
            b = ComplexA([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))
        except SingletonError as err:
            assert err.existing is a
        del a
        b = ComplexA([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))

    def test_slitting(self):
        d1, d2 = self.d1, self.d2
        a = ComplexA([d1, d1, ~d2, '+', d2, ~d1, '+', d2], list('(.(+))+.'))
        assert len(list(a.split())) == 2

if __name__ == '__main__':
    unittest.main()

