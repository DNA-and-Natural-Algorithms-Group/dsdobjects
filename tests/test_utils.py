# -*- coding: utf-8 -*-
#
# tests for dsdobjects.base_classes
#
# Written by Stefan Badelt (badelt@caltech.edu)
#

import unittest
import dsdobjects.utils as utils

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class UtilityTests(unittest.TestCase):
    def setUp(self):
        pass

    def test_make_pair_table(self):
        inp = list('(((...)))')
        res = utils.make_pair_table(inp)
        exp = [[(0, 8), (0, 7), (0, 6), None, None, None, (0, 2), (0, 1), (0, 0)]]
        rev = utils.pair_table_to_dot_bracket(exp)
        self.assertEqual(res, exp)
        self.assertEqual(inp, rev)

        inp = list('(((+...)))')
        res = utils.make_pair_table(inp)
        exp = [[(1, 5), (1, 4), (1, 3)], [None, None, None, (0, 2), (0, 1), (0, 0)]]
        rev = utils.pair_table_to_dot_bracket(exp)
        self.assertEqual(res, exp)
        self.assertEqual(inp, rev)

        res = utils.make_pair_table('((((+))).)')
        exp = [[(1, 4), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), None, (0, 0)]]
        self.assertEqual(res, exp)
        
        inp = list('((((&))).)')
        res = utils.make_pair_table(inp, strand_break='&')
        exp = [[(1, 4), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), None, (0, 0)]]
        rev = utils.pair_table_to_dot_bracket(exp, strand_break='&')
        self.assertEqual(res, exp)
        self.assertEqual(inp, rev)

        with self.assertRaises(utils.DSDUtilityError):
            res = utils.make_pair_table('((((+)).)')

        with self.assertRaises(utils.DSDUtilityError):
            res = utils.make_pair_table('((((&))).)')

    def test_make_loop_index(self):
        struct = '(((...)))'
        pt = utils.make_pair_table(struct)
        exp1 = [[1,2,3,3,3,3,3,2,1]]
        exp2 = set([0])
        out1, out2 = utils.make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

        struct = '..(((...)))..'
        pt = utils.make_pair_table(struct)
        exp1 = [[0,0,1,2,3,3,3,3,3,2,1,0,0]]
        exp2 = set([0])
        out1, out2 = utils.make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

        struct = '.((.((...))..((...).))).'
        pt = utils.make_pair_table(struct)
        exp1 = [[0,1,2,2,3,4,4,4,4,4,3,2,2,5,6,6,6,6,6,5,5,2,1,0]]
        exp2 = set([0])
        out1, out2 = utils.make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

        struct = '.((.((...)).+.((...).))).'
        pt = utils.make_pair_table(struct)
        exp1 = [[0,1,2,2,3,4,4,4,4,4,3,2],[2,5,6,6,6,6,6,5,5,2,1,0]]
        exp2 = set([0,2])
        out1, out2 = utils.make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

        struct = '.((.((...))+..((...).))).'
        pt = utils.make_pair_table(struct)
        exp1 = [[0,1,2,2,3,4,4,4,4,4,3],[2,2,5,6,6,6,6,6,5,5,2,1,0]]
        exp2 = set([0,2])
        out1, out2 = utils.make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

        struct = '.((.((...))+((...).))).'
        pt = utils.make_pair_table(struct)
        exp1 = [[0,1,2,2,3,4,4,4,4,4,3],[5,6,6,6,6,6,5,5,2,1,0]]
        exp2 = set([0,2])
        out1, out2 = utils.make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

        with self.assertRaises(utils.DSDUtilityError):
            struct = '..+..'
            pt = utils.make_pair_table(struct)
            out1, out2 = utils.make_loop_index(pt)

        with self.assertRaises(utils.DSDUtilityError):
            struct = '(.)+.(..)'
            pt = utils.make_pair_table(struct)
            out1, out2 = utils.make_loop_index(pt)

        with self.assertRaises(utils.DSDUtilityError):
            struct = '((..((.))+.(..).+((...))))'
            pt = utils.make_pair_table(struct)
            out1, out2 = utils.make_loop_index(pt)

    def test_split_complex(self):
        se = list('CCCTTTGGG')
        ss = list('(((...)))')
        ps = utils.make_lol_sequence(se)
        pt = utils.make_pair_table(ss)
        exp = [(se,ss)]
        self.assertEqual(utils.split_complex(ps, pt), exp)

        se = list('CCCT+TTGGG')
        ss = list('(((.+..)))')
        ps = utils.make_lol_sequence(se)
        pt = utils.make_pair_table(ss)
        exp = [(se,ss)]
        self.assertEqual(utils.split_complex(ps, pt), exp)

        se = list('CCCT+AAA+TTGGG')
        ss = list('(((.+...+..)))')
        ps = utils.make_lol_sequence(se)
        pt = utils.make_pair_table(ss)
        exp = [(list('CCCT+TTGGG'),list('(((.+..)))')),(list('AAA'),list('...'))]
        self.assertEqual(sorted(utils.split_complex(ps, pt)), sorted(exp))

        se = list('CCCT+AAA+TTTT+TTGGG')
        ss = list('(((.+...+....+..)))')
        ps = utils.make_lol_sequence(se)
        pt = utils.make_pair_table(ss)
        exp = [(list('CCCT+TTGGG'),list('(((.+..)))')),
               (list('AAA'),list('...')), 
               (list('TTTT'),list('....'))]
        self.assertEqual(sorted(utils.split_complex(ps, pt)), sorted(exp))

        se = list('CCCT+AAA+GCGC+TTTT+TTGGG')
        ss = list('(((.+...+)()(+....+..)))')
        ps = utils.make_lol_sequence(se)
        pt = utils.make_pair_table(ss)
        exp = [(list('CCCT+GCGC+TTGGG'),list('(((.+)()(+..)))')),
               (list('AAA'),list('...')), 
               (list('TTTT'),list('....'))]
        self.assertEqual(sorted(utils.split_complex(ps, pt)), sorted(exp))

        se = list('CCCT+AAA+GCGC+TTTT+TTGGG')
        ss = list('(((.+.(.+)()(+.)..+..)))')
        ps = utils.make_lol_sequence(se)
        pt = utils.make_pair_table(ss)
        exp = [(list('CCCT+TTGGG'),list('(((.+..)))')),
               (list('AAA+GCGC+TTTT'),list('.(.+)()(+.)..'))] 
        self.assertEqual(sorted(utils.split_complex(ps, pt)), sorted(exp))

if __name__ == '__main__':
    unittest.main()
