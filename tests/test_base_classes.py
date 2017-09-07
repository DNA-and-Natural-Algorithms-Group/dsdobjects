# -*- coding: utf-8 -*-
#
# tests for dsdobjects.base_classes
#
# Written by Stefan Badelt (badelt@caltech.edu)
#

import unittest

import dsdobjects.base_classes as bc

SKIP = False

class DummyDomain(bc.DL_Domain):
    # inherit and change the global variables ...
    # btw. doing it like this can cause troubles, if you
    # have multiple objects that get initialized like this,
    # they will all modify the bc.DL_Domain... 
    bc.DL_Domain.DTYPE_CUTOFF = 7
    bc.DL_Domain.SHORT_DOM_LEN = 6
    bc.DL_Domain.LONG_DOM_LEN = 15

    bc.DL_Domain.MEMORY = dict()

    def __init__(self, *kargs, **kwargs):
        super(DummyDomain, self).__init__(*kargs, **kwargs)

    @property
    def complement(self):
        # If we initialize the complement, we need to know the class.
        if self._complement is None:
            cname = self._name[:-1] if self.is_complement else self._name + '*'
            if cname in bc.DL_Domain.MEMORY:
                self._complement = bc.DL_Domain.MEMORY[cname]
            else :
                self._complement = DummyDomain(cname, self.dtype, self.length)
        return self._complement

@unittest.skipIf(SKIP, "skipping tests")
class Test_SequenceConstraint(unittest.TestCase):
    def setUp(self):
        pass

    def test_add_constraints(self):
        self.seq1 = 'ACGTAGTTC'
        self.seq2 = 'NNNNNNNNN'
        self.seq3 = 'RRRRRRRRR'
        self.seq4 = 'HHHHHHHHH'
        self.seq5 = 'AAAAAAAAA'
        self.seq6 = 'TGYRTYRRG'
        self.seq7 = reversed('TGYRTYRRG')

        con1 = bc.SequenceConstraint(self.seq1)
        con2 = bc.SequenceConstraint(self.seq2)
        con3 = bc.SequenceConstraint(self.seq3)
        con4 = bc.SequenceConstraint(self.seq4)
        con5 = bc.SequenceConstraint(self.seq5)
        con6 = bc.SequenceConstraint(self.seq6)
        con7 = bc.SequenceConstraint(self.seq7)

        self.assertEqual(con1, con1 + con2)

        self.assertEqual('AAAAAAAAA', str(con3 + con4))
        self.assertEqual(con5, con3 + con4)

        with self.assertRaises(bc.DSDObjectsError):
            con1 + con3

        self.assertEqual(str(con2), con2.complement)
        self.assertEqual('TGYRTYRRG', con1.complement)
        self.assertEqual(con7, ~con1)
        self.assertEqual(len(con1), 9)

@unittest.skipIf(SKIP, "skipping tests")
class DL_Domain_Test(unittest.TestCase):
    def setUp(self): 
        pass

    def tearDown(self): 
        bc.clear_memory()

    def test_initialization_memory(self):
        self.assertDictEqual(bc.DL_Domain.MEMORY, dict())
        with self.assertRaises(bc.DSDObjectsError):
            foo = bc.DL_Domain('foo')
        with self.assertRaises(bc.DSDObjectsError):
            foo = bc.DL_Domain('foo', dtype = 'blubb')
        with self.assertRaises(bc.DSDObjectsError):
            foo = bc.DL_Domain('foo', dtype = 'long', length = 5)
        self.assertDictEqual(bc.DL_Domain.MEMORY, dict())

    def test_domain_length(self):
        foo = bc.DL_Domain('foo', dtype = 'short')
        with self.assertRaises(bc.DSDDuplicationError):
            bar = bc.DL_Domain('foo', dtype = 'long')
        bar = bc.DL_Domain('bar', dtype = 'long')

        self.assertEqual(len(foo), bc.DL_Domain.SHORT_DOM_LEN)
        self.assertEqual(len(bar), bc.DL_Domain.LONG_DOM_LEN)

        foo2 = bc.DL_Domain('foo2', length = bc.DL_Domain.DTYPE_CUTOFF)
        bar2 = bc.DL_Domain('bar2', length = bc.DL_Domain.DTYPE_CUTOFF+1)

        self.assertEqual(foo2.dtype, 'short')
        self.assertEqual(bar2.dtype, 'long')

    def test_complementarity(self):
        foo = bc.DL_Domain('foo', length = 7)
        
        with self.assertRaises(bc.DSDObjectsError):
            fooC = ~foo

        fooC= bc.DL_Domain('foo*', length = 7)

        self.assertTrue(foo == ~fooC)

    def test_dummy_domain(self):
        foo2 = DummyDomain('foo2', length = 7)
        bar2 = DummyDomain('bar2', length = 8)

        self.assertEqual(foo2.dtype, 'short')
        self.assertEqual(bar2.dtype, 'long')

        self.assertTrue(foo2 != bar2)
        self.assertTrue(str(foo2), 'foo2')

    def test_dummy_complementarity(self):
        foo = DummyDomain('foo', length = 7)

        fooC = ~foo
        self.assertEqual(foo, ~fooC)
        self.assertEqual(fooC, ~foo)
        self.assertEqual(len(fooC), 7)
        self.assertEqual(fooC.dtype, 'short')

        x = ~foo
        self.assertEqual(x, fooC)
        self.assertTrue(x is fooC)

@unittest.skipIf(SKIP, "skipping tests")
class SL_Domain_Test(unittest.TestCase):
    def setUp(self): 
        self.foo = DummyDomain('foo', length = 7)
        self.bar = DummyDomain('bar', length = 10)

    def tearDown(self): 
        bc.clear_memory()

    def test_initialization(self):
        foo = bc.SL_Domain(self.foo, 'HHHHHHH')
        foo1 = bc.SL_Domain(self.foo, 'HHHHHHH', variant='1')
        foo2 = bc.SL_Domain(self.foo, 'NNNNNNN', variant='2')
        
        self.assertEqual(foo, foo1)
        self.assertFalse(foo is foo1)

        with self.assertRaises(bc.DSDObjectsError):
            bar = ~foo

        self.assertEqual(sorted(foo.variants), sorted([foo, foo1, foo2]))
        self.assertEqual(foo.length, 7)
        self.assertEqual(len(foo), 7)

        with self.assertRaises(bc.DSDObjectsError):
            foC = bc.SL_Domain(~self.foo, 'N')

        # Initialize a vaild complement domain.
        foC = bc.SL_Domain(~self.foo, 'NNNNNNN')
        self.assertEqual(foC.variants, [foC])
        # Get all complements
        self.assertEqual(sorted(~foC), sorted([foo, foo1, foo2]))
        self.assertEqual(~foo, [foC])
        self.assertEqual(~foo1, [foC])
        self.assertEqual(~foo2, [foC])

@unittest.skipIf(SKIP, "skipping tests")
class DSD_ComplexObjectTest(unittest.TestCase):
    def setUp(self):
        self.d1 = DummyDomain('d1', length=5)
        self.d2 = DummyDomain('d2', length=5)
        self.d3 = DummyDomain('d3', length=5)
        self.d1c = ~self.d1
        self.d2c = ~self.d2
        self.d3c = ~self.d3

    def tearDown(self):
        bc.clear_memory()

    def test_DSD_ComplexInit(self):
        foo = bc.DSD_Complex(sequence=list('RNNNY'), structure=list('(...)'), name='a')
        bar = bc.DSD_Complex(sequence=list('RNANY'), structure=list('(...)'), name='a')
        self.assertIsInstance(foo, bc.DSD_Complex)

        self.assertEqual(foo.sequence, list('RNNNY'))
        self.assertEqual(foo.structure, list('(...)'))
        self.assertEqual(foo.lol_sequence, [list('RNNNY')])
        self.assertEqual(foo.rotate_once(), foo)
        for r in foo.rotate():
            self.assertEqual(r.sequence, list('RNNNY'))
            self.assertEqual(r.structure, list('(...)'))

    def test_DSD_ComplexDomains(self):
        foo = bc.DSD_Complex(sequence=
                [self.d1, self.d2, self.d3, '+', self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2],
                structure=list('..(+(+))..'))

        self.assertEqual(foo.sequence, 
                [self.d1, self.d2, self.d3, '+', self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2])
        self.assertEqual(foo.lol_sequence, 
                [[self.d1, self.d2, self.d3], [self.d1], [self.d1c, self.d3c, self.d1c, self.d2]])

        bc.clear_memory()
        bar = bc.DSD_Complex(sequence=
                [self.d1c, self.d3c, self.d1c, self.d2, '+', self.d1, self.d2, self.d3, '+', self.d1], 
                structure=list('((..+..)+)'))

        self.assertEqual(foo, bar)
        self.assertTrue(foo == bar)

    def test_properties(self):
        foo = bc.DSD_Complex( sequence=[self.d1, self.d2, self.d3, '+',
            self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2], 
            structure=list('..(+(+))..'))
        with self.assertRaises(bc.DSDDuplicationError):
            bar = bc.DSD_Complex( sequence=[self.d1, self.d2, self.d3, '+',
                self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2], 
                structure=list('..(+(+))..'))

        # quick fix to rotate foo into canonical form...
        foo.rotate_once()
        tuple(map(str, foo.sequence)), tuple(foo.structure)
        canon = tuple((tuple(map(str, foo.sequence)), tuple(foo.structure)))
        foo.rotate_once()
        foo.rotate_once()

        self.assertEqual(foo.canonical_form, canon)
        with self.assertRaises(TypeError):
            foo.canonical_form[0] = 'cannot change canonical form!'

        self.assertEqual(foo.domains, [self.d1, self.d1c, self.d2, self.d3, self.d3c])
        self.assertEqual(map(str,foo.domains), ['d1', 'd1*', 'd2', 'd3', 'd3*'])

        # Make sure pair table is returned and immutable.
        self.assertEqual(foo.pair_table, 
                [[None, None, (2, 1)], [(2, 0)], [(1, 0), (0, 2), None, None]])
        pt = foo.pair_table
        self.assertTrue(foo.pair_table == pt)
        pt[0][2] = None
        self.assertFalse(foo.pair_table == pt)

        self.assertEqual(foo.get_domain((0,0)), self.d1)
        self.assertEqual(foo.get_paired_loc((0,0)), None)
        self.assertEqual(foo.get_domain((1,0)), self.d1)
        self.assertEqual(foo.get_paired_loc((1,0)), (2,0))
        self.assertEqual(foo.get_domain((2,2)), self.d1c)
        self.assertEqual(foo.get_paired_loc((2,2)), None)
        with self.assertRaises(IndexError):
            foo.get_domain((2,9))
        with self.assertRaises(IndexError):
            foo.get_paired_loc((2,9))

        self.assertEqual(foo.exterior_domains, [(0,0),(0,1),(2,2),(2,3)])

    def test_sorting(self):
        foo = bc.DSD_Complex(sequence=[self.d1, self.d2, self.d3, '+',
            self.d1, '+', self.d1c, self.d3c, self.d1c, self.d3c], structure=list('..(+(+))..'))
        bar = bc.DSD_Complex(sequence=[self.d1, self.d2, self.d3, '+',
            self.d1, '+', self.d1c, self.d3c, self.d1c, self.d3c], structure=list('..(+(+..))'))

        self.assertEqual(sorted([foo, bar]), sorted([bar, foo]))

    def test_sanitychecks(self):
        with self.assertRaises(bc.DSDObjectsError):
            # Unbalanced dot-bracket string
            foo = bc.DSD_Complex(sequence=[self.d1, self.d2, self.d3, '+', self.d1, '+', 
                self.d1c, self.d3c, self.d1c, self.d2], structure=list('..(+(+))).'))

        foo = bc.DSD_Complex( sequence=[self.d1, self.d2, self.d3, '+', self.d1, '+', 
            self.d1c, self.d3c, self.d1c, self.d2], structure=list('.((+(+))).'))

        self.assertTrue(foo.is_connected)
        with self.assertRaises(bc.DSDObjectsError):
            foo.is_domainlevel_complement

        foo = bc.DSD_Complex( sequence=[self.d1, self.d2, self.d3, '+', self.d1, '+', 
            self.d1c, self.d3c, self.d1c, self.d2], structure=list('(.(+.+.)).'))

        self.assertTrue(foo.is_domainlevel_complement)
        self.assertFalse(foo.is_connected)

    def test_names(self):
        foo = bc.DSD_Complex( sequence=[self.d1, self.d2, self.d3, '+',
            self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2], structure=list('..(+(+))..'))
        bar = bc.DSD_Complex( sequence=[self.d1, self.d2, self.d3, '+',
            self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2], structure=list('(.(+(+))).'))

        self.assertEqual(foo.name, 'cplx0')
        self.assertEqual(bar.name, 'cplx1')

        foo.name = 'foo'
        self.assertEqual(foo.name, 'foo')
        self.assertSetEqual(set(bc.DSD_Complex.NAMES.keys()), set(['foo', 'cplx1']))

        foo.name = 'bar'
        self.assertEqual(foo.name, 'bar')
        self.assertSetEqual(set(bc.DSD_Complex.NAMES.keys()), set(['bar', 'cplx1']))

        with self.assertRaises(bc.DSDObjectsError):
            bar.name = 'bar'

    def test_rotations(self):
        foo = bc.DSD_Complex(sequence=[self.d1, self.d2, self.d3, '+',
            self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2],
                structure=list('..(+(+))..')) 
        self.assertEqual(foo.rotate_once(), foo)
        self.assertEqual(foo.sequence, 
                [self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2, '+', 
                    self.d1, self.d2, self.d3])
        self.assertEqual(foo.structure, list('(+)(..+..)'))

        for r in foo.rotate():
            self.assertEqual(r, foo)
        self.assertEqual(foo.sequence, [self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2, '+', 
                             self.d1, self.d2, self.d3])
        self.assertEqual(foo.structure, list('(+)(..+..)'))

@unittest.skipIf(SKIP, "skipping tests")
class DSD_ReactionTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        bc.clear_memory()

    def test_initialize(self):
        with self.assertRaises(bc.DSDObjectsError):
            x = bc.DSD_Reaction(['A','B','C'],['D','A','E'], rtype = 'branch-3way')

        A = bc.DSD_Complex(list('NNNNNNN'), list('((...))'), name='A')
        B = bc.DSD_Complex(list('NNNNNNN'), list('.(...).'), name='B')

        x = bc.DSD_Reaction([A, B],[B, B], rtype = 'branch-3way')
        with self.assertRaises(bc.DSDDuplicationError):
            y = bc.DSD_Reaction([A, B],[B, B], rtype = 'branch-3way')

        bc.clear_memory()
        y = bc.DSD_Reaction([A, B],[B, B], rtype = 'branch-3way', rate=.5)
        z = bc.DSD_Reaction([A],[B], rtype='bind11')

        self.assertEqual('A + B -> B + B', str(x))
        A.name = 'Z'
        self.assertEqual('Z + B -> B + B', str(x))
        self.assertEqual('Z + B -> B + B', str(y))

        self.assertEqual(x, y)

        with self.assertRaises(bc.DSDDuplicationError):
            y = bc.DSD_Reaction([A, B],[B, B], rtype = 'branch-3way')

        self.assertEqual(x, y)

        self.assertTrue(x != z)
        self.assertEqual(y.rate, .5)
        self.assertEqual(y.rateunits, '/M/s')


if __name__ == '__main__':
    unittest.main()

