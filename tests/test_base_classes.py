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
from dsdobjects.base_classes import DomainS, StrandS, ComplexS

SKIP = False
SKIP_TOXIC = True

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

    def test_initialization_03(self):
        x = DomainS('a', 15)
        y = DomainS('a*')
        assert len(x) == len(y)
        assert x is ~y
        assert y is ~x

        x = DomainS('b*', 5)
        with self.assertRaises(SingletonError):
            y = DomainS('b')
        with self.assertRaises(SingletonError):
            y = DomainS('b', 10)
        y = DomainS('b', 5)
        assert len(x) == len(y)
        assert x is ~y
        assert y is ~x

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

    def test_dl_initialization_01(self):
        d1, d2, d3 = self.d1, self.d2, self.d3
        d1c, d2c, d3c = self.d1c, self.d2c, self.d3c

        foo = ComplexS(sequence =  [d1,  d2,  d3,  '+', d1,  '+', d1c, d3c, d1c, d2],
                       structure = ['.', '.', '(', '+', '(', '+', ')', ')', '.', '.'],
                       name = 'foo')
        try:
            bar = ComplexS(sequence = [d1c, d3c, d1c, d2, '+', d1, d2, d3, '+', d1], 
                       structure = list('((..+..)+)'))
        except SingletonError as err:
            bar = err.existing
        self.assertTrue(foo is bar)
        self.assertTrue(foo == bar)

        seq = list(foo.sequence)
        assert seq[3] == '+'
        seq[3] = '-'
        assert seq != list(foo.sequence)

        with self.assertRaises(SingletonError):
            bar = ComplexS(sequence = [d1c, d3c, d1c, d2, '+', d1, d2, d3, '+', d1], 
                           structure = list('((..+..)+)'), 
                           name = 'bar')

        assert foo.canonical_form == (('d1', '+', 'd1*', 'd3*', 'd1*', 'd2', '+', 'd1', 'd2', 'd3'), 
                                      ('(', '+', ')', '(', '.', '.', '+', '.', '.', ')'))
        assert foo.turns == 2 
        foo.turns = 0 # rotate foo into canonical form ...
        assert foo.kernel_string == 'd1( + ) d3*( d1* d2 + d1 d2 )'
        foo.turns += 1
        assert foo.kernel_string == 'd1*( d3*( d1* d2 + d1 d2 ) + )'
        foo.turns += 1
        assert foo.kernel_string == 'd1 d2 d3( + d1( + ) ) d1* d2'
        foo.turns += 1
        assert foo.kernel_string == 'd1( + ) d3*( d1* d2 + d1 d2 )'
        assert foo.turns == 0

    def test_rotations(self):
        d1, d2, d3 = self.d1, self.d2, self.d3
        d1c, d2c, d3c = self.d1c, self.d2c, self.d3c

        foo = ComplexS(sequence =  [d1,  d2,  d3,  '+', d1,  '+', d1c, d3c, d1c, d2],
                       structure = ['.', '.', '(', '+', '(', '+', ')', ')', '.', '.'],
                       name = 'foo')
        loc0 = (0, 2)
        loc1 = (1, 0)
        loc2 = (2, 3)

        ploc0 = foo.get_paired_loc(loc0)
        dom0 = foo.get_domain(loc0)
        assert ploc0 == (2, 1)
        assert dom0 is d3
        ploc1 = foo.get_paired_loc(loc1)
        dom1 = foo.get_domain(loc1)
        assert dom1 is d1
        assert ploc1 == (2, 0)
        ploc2 = foo.get_paired_loc(loc2)
        dom2 = foo.get_domain(loc2)
        assert ploc2 == None
        assert dom2 is d2

        for e, (x, y) in enumerate(foo.rotate_pt()):
            l0 = foo.rotate_pairtable_loc(loc0, e)
            l1 = foo.rotate_pairtable_loc(loc1, e)
            l2 = foo.rotate_pairtable_loc(loc2, e)
            assert x[l0[0]][l0[1]] is dom0
            assert x[l1[0]][l1[1]] is dom1
            assert x[l2[0]][l2[1]] is dom2

 
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
        self.assertEqual(list(foo.pair_table), pt)
        pt = list(foo.pair_table)
        pt[1][0] = None
        self.assertFalse(list(foo.pair_table) == pt)

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
        assert str(a) == ComplexS.PREFIX + '1'
        b = ComplexS(list('A'), list('.'))
        assert b.kernel_string == 'A'
        assert str(b) == ComplexS.PREFIX + '2'
        assert a is ComplexS(None, None, name = ComplexS.PREFIX + '1')
    
    def test_initialization_02(self):
        d1, d2 = self.d1, self.d2
        a = ComplexS([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))
        assert a is ComplexS([d2, ~d1, '+', d1, d1, ~d2], list('((+).)'), name = a.name)
        assert a is ComplexS(None, None, name = a.name)

    def test_initialization_03(self):
        ComplexS.ID = 1
        ComplexS.PREFIX = 'c'
        c1 = ComplexS(list('ABCDEFG'), list('.......'), name = 'c1') 
        with self.assertRaises(SingletonError):
            c2 = ComplexS(list('BACEFGD'), list('.......'))
    
    def test_exceptions(self):
        d1, d2 = self.d1, self.d2
        a = ComplexS([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))
        try:
            b = ComplexS([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))
        except SingletonError as err:
            assert err.existing is a
        del a
        b = ComplexS([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'))

    def test_splitting(self):
        d1, d2 = self.d1, self.d2
        a = ComplexS([d1, d1, ~d2, '+', d2, ~d1, '+', d2], list('(.(+))+.'))
        assert len(list(a.split())) == 2

    def test_split_exception(self):
        d1, d2 = self.d1, self.d2
        ComplexS.ID = 1
        ComplexS.PREFIX = 'c'
        c1 = ComplexS([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'), name = 'c1')
        c2 = ComplexS([d2], list('.'), name = 'c2')
        a = ComplexS([d1, d1, ~d2, '+', d2, ~d1, '+', d2], list('(.(+))+.'), name = 'a')
        with self.assertRaises(SingletonError):
            list(a.split())

    def test_common_practice(self):
        d1, d2 = self.d1, self.d2
        ComplexS.ID = 1
        ComplexS.PREFIX = 'c'
        c1 = ComplexS([d1, d1, ~d2, '+', d2, ~d1], list('(.(+))'), name = 'c1')
        try:
            c2 = ComplexS([d2], list('.'))
        except SingletonError as err:
            assert err.existing is None

@unittest.skipIf(SKIP, "skipping tests.")
class TestStrandS(unittest.TestCase):
    def setUp(self):
        self.d1 = DomainS('d1', length = 20)
        self.d2 = DomainS('d2', length = 5)
        self.d3 = DomainS('d3', length = 15)

    def tearDown(self):
        clear_singletons(DomainS)
        clear_singletons(StrandS)

    def test_strand_initialization_01(self):
        d1, d2, d3 = self.d1, self.d2, self.d3
        foo = StrandS(sequence = [d1, d2, ~d3], name = 'foo')
        assert str(foo) == 'foo'
        assert foo.size == 1
        assert foo.canonical_form == (('d1', 'd2', 'd3*'), tuple('***'))
        assert list(foo.sequence) == [d1, d2, ~d3]
        assert foo.structure is None


if __name__ == '__main__':
    unittest.main()

