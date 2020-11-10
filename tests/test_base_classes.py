#
# tests/test_base_classes.py
#   - copy and/or modify together with dsdobjects/base_classes.py
#
import logging
logger = logging.getLogger('dsdobjects')
logger.setLevel(logging.INFO)
import unittest

from dsdobjects import (SecondaryStructureError,
                        SingletonError,
                        show_singletons,
                        clear_singletons)
from dsdobjects.base_classes import DomainS, ComplexS

SKIP = False
SKIP_TOXIC = True
SKIP_DEPRECATED = True

class MyDomainS(DomainS):
    pass

@unittest.skipIf(SKIP, "skipping tests.")
class TestSingletonInheritance(unittest.TestCase):
    def tearDown(self):
        clear_singletons(MyDomainS)
        clear_singletons(DomainS)

    @unittest.skipIf(SKIP_TOXIC, "skipping toxic tests.")
    def test_garbage_collection(self):
        # This test must fail to see if garbage collection works.
        # Spoiler: it doesn't that is why we need tearDown()
        a = DomainS('a', 10)
        a = MyDomainS('a', 10)
        assert False

    def test_initialization_01(self):
        a = DomainS('a', 10)
        b = MyDomainS('a', 10)
        assert a is not b
        assert a == b

    def test_initialization_02(self):
        assert DomainS('a', 10) == ~MyDomainS('a*', 10)
        assert DomainS('a', 10) != ~MyDomainS('a*', 15)

@unittest.skipIf(SKIP, "skipping tests.")
class TestSingletonDomain(unittest.TestCase):
    def tearDown(self):
        clear_singletons(DomainS)

    @unittest.skipIf(SKIP_TOXIC, "skipping toxic tests.")
    def test_garbage_collection(self):
        # This test must fail to see if garbage collection works.
        # Spoiler: it doesn't that is why we need tearDown()
        a = DomainS('a', 10)
        assert False

    def test_initialization_01(self):
        # Not enough parameters to initialize
        with self.assertRaises(SingletonError):
            a = DomainS('a')

        a = DomainS('a', 10)
        assert a is DomainS('a')

        b = ~a
        del a, b # Delete a, b and its reference!

        # Check complementarity on anonymos objects.
        assert DomainS('a', 15) is ~DomainS('a*', 15)

        a = DomainS('a', 10)
        with self.assertRaises(SingletonError):
            b = DomainS('a', 12)

    def test_initialization_02(self):
        a = DomainS('a', 10)
        b = DomainS('b', 10)
        assert a is DomainS('a')
        del a
        ## Check complementarity on anonymos objects.
        assert DomainS('a', 10) is ~DomainS('a*', 10)

    def test_immutable_forms(self):
        a = DomainS('a', 15)
        b = DomainS('b', 15)

        # Forbidden change of immutable attributes.
        with self.assertRaises(SingletonError):
            a.length = 10
        assert len(a) == 15

        with self.assertRaises(SingletonError):
            a.name = 'b'
        assert a.name == 'a'

        with self.assertRaises(SingletonError):
            a = DomainS('b', 10)

@unittest.skipIf(SKIP, "skipping tests.")
class TestAutomaticDomain(unittest.TestCase):
    def tearDown(self):
        clear_singletons(DomainS)

    @unittest.skipIf(SKIP_TOXIC, "skipping toxic tests.")
    def test_garbage_collection(self):
        # This test must fail to see if garbage collection works.
        # Spoiler: it doesn't that is why we need tearDown()
        a = DomainS('a', 10)
        assert False

    def test_initialization_01(self):
        # Not enough parameters to initialize
        with self.assertRaises(SingletonError):
            a = DomainS('a')
        a = DomainS('a', 10)
        assert a is DomainS('a')
        assert DomainS('a', 10) is ~DomainS('a*', 10)

    def test_initialization_02(self):
        DomainS.ID = 1
        a = DomainS(dtype = 'short')
        b = DomainS(prefix = 'a', length = 10)
        c = DomainS(prefix = 'a', length = 10)
        assert (a.name, a.length) == ('d1', 5)
        assert (b.name, b.length) == ('a2', 10)
        assert (c.name, c.length) == ('a3', 10)

    def test_immutable_forms(self):
        a = DomainS('a', 15)
        b = DomainS('b', 15)

        # Forbidden change of immutable attributes.
        with self.assertRaises(SingletonError):
            a.length = 10
        assert len(a) == 15

        with self.assertRaises(SingletonError):
            a.name = 'b'
        assert a.name == 'a'

        with self.assertRaises(SingletonError):
            a = DomainS('b', 10)

@unittest.skipIf(SKIP, "skipping tests.")
class TestSingletonComplex(unittest.TestCase):
    def setUp(self):
        self.d1 = DomainS('d1', 5)
        self.d2 = DomainS('d2', 5)
        self.d3 = DomainS('d3', 5)
        self.d1c = ~self.d1
        self.d2c = ~self.d2
        self.d3c = ~self.d3

    def tearDown(self):
        clear_singletons(DomainS)
        clear_singletons(ComplexS)

    def test_nl_initialization_01(self):
        foo = ComplexS(sequence = list('RNNNY'), 
                       structure = list('(...)'), 
                       name = 'a')
        bar = ComplexS(sequence = list('RNNNY'), structure = list('(...)'), name = 'a')
        assert foo == bar

        with self.assertRaises(SingletonError):
            bar = ComplexS(sequence = list('RNNNY'), structure = list('(...)'), name = 'b')

        with self.assertRaises(SingletonError):
            bar = ComplexS(sequence = list('RNANY'), structure = list('(...)'), name = 'a')

        self.assertEqual(foo.sequence, list('RNNNY'))
        self.assertEqual(foo.structure, list('(...)'))
        self.assertEqual(foo.strand_table, [list('RNNNY')])
        self.assertEqual(next(foo.rotate(1)), (foo.sequence, foo.structure))

        del foo, bar
        bar = ComplexS(sequence = list('RNANY'), structure = list('(...)'), name = 'a')

    def test_dl_initialization_01(self):
        d1, d2, d3 = self.d1, self.d2, self.d3
        d1c, d2c, d3c = self.d1c, self.d2c, self.d3c

        foo = ComplexS(sequence =  [d1,  d2,  d3,  '+', d1,  '+', d1c, d3c, d1c, d2],
                       structure = ['.', '.', '(', '+', '(', '+', ')', ')', '.', '.'],
                       name = 'fred')
        try:
            bar = ComplexS(sequence = [d1c, d3c, d1c, d2, '+', d1, d2, d3, '+', d1], 
                       structure = list('((..+..)+)'))
        except SingletonError as err:
            bar = err.existing

        assert foo is bar
        self.assertEqual(foo, bar)
        self.assertTrue(foo == bar)

        with self.assertRaises(SingletonError):
            bar = ComplexS(sequence = [d1c, d3c, d1c, d2, '+', d1, d2, d3, '+', d1], 
                           structure = list('((..+..)+)'), 
                           name = 'barney')

        #for e, (x, y) in enumerate(bar.rotate()):
        #    print(e, list(map(str, x)), y)

        assert foo.canonical_form == (('d1', '+', 'd1*', 'd3*', 'd1*', 'd2', '+', 'd1', 'd2', 'd3'), ('(', '+', ')', '(', '.', '.', '+', '.', '.', ')'))
        assert foo.turns == 1 

        foo.turns = 0 # rotate foo into canonical form ...
        assert foo.kernel_string == 'd1( + ) d3*( d1* d2 + d1 d2 )'
        foo.turns += 1
        assert foo.kernel_string == 'd1 d2 d3( + d1( + ) ) d1* d2'
        foo.turns += 1
        assert foo.kernel_string == 'd1*( d3*( d1* d2 + d1 d2 ) + )'
        foo.turns += 1
        assert foo.kernel_string == 'd1( + ) d3*( d1* d2 + d1 d2 )'
        assert foo.turns == 0

    def test_init_disconnected(self):
        d1, d2, d3 = self.d1, self.d2, self.d3
        d1c, d2c, d3c = self.d1c, self.d2c, self.d3c

        foo = ComplexS(sequence =  [d1,  d2,  d3,  '+', d1,  '+', d1c, d3c, d1c, d2],
                       structure = list('...+.+....'), name = 'foo')
        assert foo.kernel_string == 'd1 d2 d3 + d1 + d1* d3* d1* d2'
        with self.assertRaises(SecondaryStructureError):
            foo.exterior_domains
        assert not foo.is_connected
        f1 = ComplexS([d1,  d2,  d3], ['.', '.', '.'], 'a')
        f2 = ComplexS([d1], ['.'], 'b')
        f3 = ComplexS([d1c, d3c, d1c, d2], ['.', '.', '.', '.'], 'c')
        for c in foo.split():
            assert c in [f1, f2, f3]

    def test_domain_properties(self):
        d1, d2, d3 = self.d1, self.d2, self.d3
        d1c, d2c, d3c = self.d1c, self.d2c, self.d3c

        foo = ComplexS(sequence = [d1, d2, d3], structure=list('...'), name = 'foo')
        self.assertEqual(foo.domains, set([d1, d2, d3]))
        del foo

        foo = ComplexS(sequence = [d1, d2, d3, '+', d1, '+', d1c, d3c, d1c, d2], 
                       structure=list('..(+(+))..'), name = 'foo')
        self.assertEqual(foo.domains, set([d1, d1c, d2, d3, d3c]))

        self.assertEqual(foo.get_domain((0,0)), d1)
        self.assertEqual(foo.get_paired_loc((0,0)), None)
        self.assertEqual(foo.get_domain((1,0)), d1)
        self.assertEqual(foo.get_paired_loc((1,0)), (2,0))
        self.assertEqual(foo.get_domain((2,2)), d1c)
        self.assertEqual(foo.get_paired_loc((2,2)), None)
        with self.assertRaises(IndexError):
            foo.get_domain((2,9))
        with self.assertRaises(IndexError):
            foo.get_paired_loc((2,9))
        with self.assertRaises(IndexError):
            foo.get_paired_loc((1,-1))
        self.assertEqual(foo.exterior_domains, [(0,0), (0,1), (2,2), (2,3)])
        self.assertEqual(foo.enclosed_domains, [])

    def test_other_properties(self):
        d1, d2, d3 = self.d1, self.d2, self.d3
        d1c, d2c, d3c = self.d1c, self.d2c, self.d3c
        foo = ComplexS(sequence = [d1, d2, d3, '+', d1, '+', d1c, d3c, d1c, d2], 
                       structure = list('..(+(+))..'), name = 'foo')
        pt = [[None, None, (2, 1)], [(2, 0)], [(1, 0), (0, 2), None, None]]
        self.assertEqual(foo.pair_table, pt)
        pt = foo.pair_table
        pt[0][2] = None
        self.assertFalse(foo.pair_table == pt)

    def test_strand_initialization_01(self):
        foo = ComplexS(sequence = list('RNNNY'), structure = None, name = 'foo')
        assert str(foo) == 'foo'
        assert foo.size == 1
        assert foo.canonical_form == (('R', 'N', 'N', 'N', 'Y'), tuple('*****'))
        assert foo.sequence == ['R', 'N', 'N', 'N', 'Y']
        assert foo.structure is None

    def test_strand_initialization_02(self):
        d1, d2, d3 = self.d1, self.d2, self.d3
        foo = ComplexS(sequence = [d1, d2, ~d3], structure = None, name = 'foo')
        assert str(foo) == 'foo'
        assert foo.size == 1
        assert foo.canonical_form == (('d1', 'd2', 'd3*'), tuple('***'))
        assert foo.sequence == [d1, d2, ~d3]
        assert foo.structure is None

@unittest.skipIf(SKIP, "skipping tests.")
class TestAutomaticComplex(unittest.TestCase):
    def setUp(self):
        self.d1 = DomainS(length = 20)
        self.d2 = DomainS(length = 5)

    def tearDown(self):
        clear_singletons(DomainS)
        clear_singletons(ComplexS)

    @unittest.skipIf(SKIP_TOXIC, "skipping toxic tests.")
    def test_garbage_collection(self):
        # This test must fail to see if garbage collection works.
        # Spoiler: it doesn't that is why we need tearDown()
        a = ComplexS(list('ABCDEFG'), list('.......')) 
        assert False

    def test_initialization_01(self):
        # Not enough parameters to initialize
        ComplexS.ID = 1
        with self.assertRaises(SingletonError):
            a = ComplexS(None, None, name = 'hello')
        a = ComplexS(list('ABCDEFG'), list('.......')) 
        assert a.kernel_string == 'A B C D E F G'
        assert str(a) == 'cplx1'
        b = ComplexS(list('A'), list('.'))
        assert b.kernel_string == 'A'
        assert str(b) == 'cplx2'
        assert a is ComplexS(None, None, name = 'cplx1')
    
    def test_initialization_02(self):
        d1, d2 = self.d1, self.d2
        a = ComplexS([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))
        assert a is ComplexS([d2, ~d1, '+', d1, d1, ~d2], list('((+).)'), name = a.name)
        assert a is ComplexS(None, None, name = a.name)

    def test_exceptions(self):
        d1, d2 = self.d1, self.d2
        a = ComplexS([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))
        try:
            b = ComplexS([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))
        except SingletonError as err:
            assert err.existing is a
        del a
        b = ComplexS([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))

    def test_slitting(self):
        d1, d2 = self.d1, self.d2
        a = ComplexS([d1, d1, ~d2, '+', d2, ~d1, '+', d2], list('(.(+))+.'))
        assert len(list(a.split())) == 2

if __name__ == '__main__':
    unittest.main()

