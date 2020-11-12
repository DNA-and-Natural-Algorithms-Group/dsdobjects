#
# tests/test_complex_utils.py
#   - copy and/or modify together with dsdobjects/complex_utils.py
#
import logging
logger = logging.getLogger('dsdobjects')
logger.setLevel(logging.INFO)
import unittest

from dsdobjects.complex_utils import (SecondaryStructureError,
                                      make_pair_table,
                                      pair_table_to_dot_bracket, 
                                      make_strand_table,
                                      strand_table_to_sequence,
                                      split_complex_db,
                                      split_complex_pt,
                                      rotate_complex_db,
                                      rotate_complex_pt,
                                      rotate_complex_once,
                                      make_loop_index)

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class TestTables(unittest.TestCase):
    def test_pair_tables(self):
        inp = '(((...)))'
        res = make_pair_table(inp)
        exp = [[(0, 8), (0, 7), (0, 6), None, None, None, (0, 2), (0, 1), (0, 0)]]
        rev = pair_table_to_dot_bracket(exp, join = True)
        self.assertEqual(res, exp)
        self.assertEqual(inp, rev)

        inp = '(((+...)))'
        res = make_pair_table(inp)
        exp = [[(1, 5), (1, 4), (1, 3)], [None, None, None, (0, 2), (0, 1), (0, 0)]]
        rev = pair_table_to_dot_bracket(exp, join = True)
        self.assertEqual(res, exp)
        self.assertEqual(inp, rev)

        res = make_pair_table('((((+))).)')
        exp = [[(1, 4), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), None, (0, 0)]]
        self.assertEqual(res, exp)
        
        inp = '((((&))).)'
        res = make_pair_table(inp, strand_break='&')
        exp = [[(1, 4), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), None, (0, 0)]]
        rev = pair_table_to_dot_bracket(exp, strand_break='&', join = True)
        self.assertEqual(res, exp)
        self.assertEqual(inp, rev)

        with self.assertRaises(SecondaryStructureError):
            res = make_pair_table('((((+)).)')

        with self.assertRaises(SecondaryStructureError):
            res = make_pair_table('((((&))).)')

    def test_strand_tables(self):
        se = 'CCCTTTGGG'
        assert make_strand_table(se) == [list(se)]
        assert strand_table_to_sequence(make_strand_table(se), join = True) == se

        se = 'CCCT+TTGGG'
        assert make_strand_table(se) == [list('CCCT'), list('TTGGG')]
        assert strand_table_to_sequence(make_strand_table(se), join = True) == se

        se = 'CCCT+TTGGG'
        assert make_strand_table(se) == [list('CCCT'), list('TTGGG')]
        assert strand_table_to_sequence(make_strand_table(se), join = True) == se

        se = 'CCCT+AAA+TTGGG'
        assert make_strand_table(se) == [list('CCCT'), list('AAA'), list('TTGGG')]
        assert strand_table_to_sequence(make_strand_table(se), join = True) == se

        se = list('CCCT+AAA+TTGGG')
        assert make_strand_table(se) == [list('CCCT'), list('AAA'), list('TTGGG')]
        assert strand_table_to_sequence(make_strand_table(se), join = False) == se

class TestLoopIndex(unittest.TestCase):
    def test_make_loop_index_00(self):
        struct = '.'
        pt = make_pair_table(struct)
        exp1 = [[0]]
        exp2 = set([0])
        out1, out2 = make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)


        struct = '(((...)))'
        pt = make_pair_table(struct)
        exp1 = [[1,2,3,3,3,3,3,2,1]]
        exp2 = set([0])
        out1, out2 = make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

        struct = '..(((...)))..'
        pt = make_pair_table(struct)
        exp1 = [[0,0,1,2,3,3,3,3,3,2,1,0,0]]
        exp2 = set([0])
        out1, out2 = make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

        struct = '.((.((...))..((...).))).'
        pt = make_pair_table(struct)
        exp1 = [[0,1,2,2,3,4,4,4,4,4,3,2,2,5,6,6,6,6,6,5,5,2,1,0]]
        exp2 = set([0])
        out1, out2 = make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

    def test_make_loop_index_01(self):
        struct = '.((.((...)).+.((...).))).'
        pt = make_pair_table(struct)
        exp1 = [[0,1,2,2,3,4,4,4,4,4,3,2],[2,5,6,6,6,6,6,5,5,2,1,0]]
        exp2 = set([0,2])
        out1, out2 = make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

        struct = '.((.((...))+..((...).))).'
        pt = make_pair_table(struct)
        exp1 = [[0,1,2,2,3,4,4,4,4,4,3],[2,2,5,6,6,6,6,6,5,5,2,1,0]]
        exp2 = set([0,2])
        out1, out2 = make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

        struct = '.((.((...))+((...).))).'
        pt = make_pair_table(struct)
        exp1 = [[0,1,2,2,3,4,4,4,4,4,3],[5,6,6,6,6,6,5,5,2,1,0]]
        exp2 = set([0,2])
        out1, out2 = make_loop_index(pt)
        self.assertEqual(out1, exp1)
        self.assertSetEqual(out2, exp2)

    def test_make_loop_index_02(self):
        with self.assertRaises(SecondaryStructureError):
            struct = '..+..'
            pt = make_pair_table(struct)
            out1, out2 = make_loop_index(pt)

        with self.assertRaises(SecondaryStructureError):
            struct = '(.)+.(..)'
            pt = make_pair_table(struct)
            out1, out2 = make_loop_index(pt)

        with self.assertRaises(SecondaryStructureError):
            struct = '((..((.))+.(..).+((...))))'
            pt = make_pair_table(struct)
            out1, out2 = make_loop_index(pt)

    def test_make_loop_index_03(self):
        struct = '..+..'
        pt = make_pair_table(struct)
        out = make_loop_index(pt, components = True)
        assert out == ([[0, 0], [0, 0]], 
                       [[0, 0], [0, 0]])

        struct = '(.)+.(..)'
        pt = make_pair_table(struct)
        out = make_loop_index(pt, components = True)
        assert out == ([[1, 1, 1], [0, 2, 2, 2, 2]], 
                       [[0, 0], [0, 0]])

        struct = '(.)(+(+(+(.)+)+)+).(..)'
        pt = make_pair_table(struct)
        out = make_loop_index(pt, components = True)
        assert out == ([[1, 1, 1, 2], [3], [4], [5, 5, 5], [4], [3], [2, 0, 6, 6, 6, 6]], 
                       [[0, 2], [2, 3], [3, 4], [4, 4], [4, 3], [3, 2], [2, 0]])

        struct = '(.)(+(+(+(+))+(.)+)+).(..)'
        pt = make_pair_table(struct)
        out = make_loop_index(pt, components = True)
        assert out == ([[1, 1, 1, 2], [3], [4], [5], [5, 4], [6, 6, 6], [3], [2, 0, 7, 7, 7, 7]], 
                       [[0, 2], [2, 3], [3, 4], [4, 5], [5, 3], [3, 3], [3, 2], [2, 0]])

        struct = '((..((.))+.(..).+((...))))'
        pt = make_pair_table(struct)
        out = make_loop_index(pt, components = True)
        assert out == ([[1, 2, 2, 2, 3, 4, 4, 4, 3], [2, 5, 5, 5, 5, 2], 
                        [6, 7, 7, 7, 7, 7, 6, 2, 1]], 
                       [[0, 2], [2, 2], [2, 0]])

        struct = '((..((.))+.(..).+(+(...)+)))'
        pt = make_pair_table(struct)
        out = make_loop_index(pt, components = True)
        assert out == ([[1, 2, 2, 2, 3, 4, 4, 4, 3], [2, 5, 5, 5, 5, 2], 
                        [6], [7, 7, 7, 7, 7], [6, 2, 1]], 
                       [[0, 2], [2, 2], [2, 6], [6, 6], [6, 0]])

        struct = '((..((.))+(+.(..)+).+(+(...)+)))'
        pt = make_pair_table(struct)
        out = make_loop_index(pt, components = True)
        assert out == ([[1, 2, 2, 2, 3, 4, 4, 4, 3], [5], [5, 6, 6, 6, 6],
                        [5, 2], [7], [8, 8, 8, 8, 8], [7, 2, 1]],
                       [[0, 2], [2, 5], [5, 5], [5, 2], [2, 7], [7, 7], [7, 0]])

class TestComplexOperations(unittest.TestCase):
    def test_split_complex_00(self):
        se = list('CCCTTTGGG')
        ss = list('(((...)))')
        ps = make_strand_table(se)
        pt = make_pair_table(ss)
        out = list(split_complex_pt(ps, pt))
        assert len(out) == 1
        c1 = ([['C', 'C', 'C', 'T', 'T', 'T', 'G', 'G', 'G']], 
              [[(0, 8), (0, 7), (0, 6), None, None, None, (0, 2), (0, 1), (0, 0)]])
        assert c1 in out

        outdb = list(split_complex_db(se, ss))
        assert (se, ss) in outdb

        se = list('CCCT+TTGGG')
        ss = list('(((.+..)))')
        ps = make_strand_table(se)
        pt = make_pair_table(ss)
        out = list(split_complex_pt(ps, pt))
        assert len(out) == 1
        outdb = list(split_complex_db(se, ss))
        assert (se, ss) in outdb

    def test_split_complex_01(self):
        se = 'CCCT+AAA+TTGGG'
        ss = '(((.+...+..)))'
        ps = make_strand_table(se)
        pt = make_pair_table(ss)
        out = list(split_complex_pt(ps, pt))
        assert len(out) == 2
        outdb = list(split_complex_db(se, ss, join = True))
        assert ('CCCT+TTGGG', '(((.+..)))') in outdb
        assert ('AAA', '...') in outdb

        se = 'CCCT+AAA+TTTT+TTGGG'
        ss = '(((.+...+....+..)))'
        ps = make_strand_table(se)
        pt = make_pair_table(ss)
        out = list(split_complex_pt(ps, pt))
        assert len(out) == 3
        outdb = list(split_complex_db(se, ss, join = True))
        assert ('CCCT+TTGGG', '(((.+..)))') in outdb
        assert ('AAA', '...') in outdb
        assert ('TTTT', '....') in outdb

    def test_split_complex_02(self):
        se = 'CCCT+AAA+GCGC+TTTT+TTGGG'
        ss = '(((.+.(.+)()(+.)..+..)))'
        ps = make_strand_table(se)
        pt = make_pair_table(ss)
        bs, bt = ps[:], pt[:]
        out = list(split_complex_pt(ps, pt))
        assert len(out) == 2
        assert bs == ps
        assert bt == pt
        outdb = list(split_complex_db(se, ss, join = True))
        assert len(outdb) == 2
        assert ('CCCT+TTGGG', '(((.+..)))') in outdb
        assert ('AAA+GCGC+TTTT', '.(.+)()(+.)..') in outdb

    def test_split_complex_03(self):
        se = 'CCCT+TC+GT+TG+C+GG'
        ss = '(((.+.(+).+.)+.+))'
        ps = make_strand_table(se)
        pt = make_pair_table(ss)
        out = list(split_complex_pt(ps, pt))
        assert len(out) == 3
        c1 = ([['T', 'C'], ['G', 'T']], 
              [[None, (1, 0)], [(0, 1), None]])
        c2 = ([['C']], 
              [[None]])
        c3 = ([['C', 'C', 'C', 'T'], ['T', 'G'], ['G', 'G']],
              [[(2, 1), (2, 0), (1, 1), None], [None, (0, 2)], [(0, 1), (0, 0)]])
        assert c1 in out
        assert c2 in out
        assert c3 in out

        outdb = list(split_complex_db(se, ss, join = True))
        assert ('TC+GT', '.(+).') in outdb
        assert ('C', '.') in outdb
        assert ('CCCT+TG+GG', '(((.+.)+))') in outdb

    def test_rotate_complex_00(self):
        se = list('CCCTTTGGG')
        ss = list('(((...)))')
        out = list(rotate_complex_db(se, ss))
        assert len(out) == 1
        out = list(rotate_complex_db(se, ss, turns = 1))
        assert len(out) == 1
        for nse, nss in out:
            assert nse == se
            assert nss == ss
        out = list(rotate_complex_db(se, ss, turns = 2))
        assert len(out) == 2

    def test_rotate_complex_01(self):
        se = list('CCCT+TTGGG')
        ss = list('(((.+..)))')
        out = list(rotate_complex_db(se, ss))
        assert len(out) == 2
        assert (se, ss) == out[0]
        out = list(rotate_complex_db(se, ss, turns = 1))
        assert len(out) == 1
        assert (se, ss) != out[0]
        out = list(rotate_complex_db(se, ss, turns = 10))
        assert len(out) == 10
        assert (se, ss) != out[0]

    def test_rotate_complex_01(self):
        se = list('CCCT+TTGGG')
        ss = list('(((.+..)))')
        out = list(rotate_complex_db(se, ss))
        assert len(out) == 2
        assert (se, ss) == out[0]
        out = list(rotate_complex_db(se, ss, turns = 1))
        assert len(out) == 1
        assert (se, ss) != out[0]
        out = list(rotate_complex_db(se, ss, turns = 10))
        assert len(out) == 10
        assert (se, ss) != out[0]

    def test_rotate_complex_once(self):
        se = list('CCCT+TTGGG')
        ss = list('(((.+..)))')
        re = list('TTGGG+CCCT')
        rs = list('..(((+))).')
        out = rotate_complex_once(se, ss)
        assert out == (re, rs)

if __name__ == '__main__':
    unittest.main()
