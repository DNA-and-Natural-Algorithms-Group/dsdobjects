# -*- coding: utf-8 -*-
#
# tests for dsdobjects.parser.pil_kernel_format
#
# Written by Stefan Badelt (badelt@caltech.edu)
#

import unittest
from pyparsing import ParseException

from dsdobjects.parser import parse_kernel_file, parse_kernel_string

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class TestKERNELparser(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_new_kernel_format(self):
        # Names starting with digits are forbidden,
        # this can easily be changed in the pil_parser.
        with self.assertRaises(ParseException):
            parse_kernel_string("length 1 = 6 ")
        with self.assertRaises(ParseException):
            parse_kernel_string("length 1ta = 6 ")

    def test_self_loops(self):
        self.assertEqual(parse_kernel_string("cplx = a( b( c ) )"),
                         [['complex', 'cplx', ['a', ['b', ['c']]]]])
        self.assertEqual(parse_kernel_string("cplx = a( b( c( ) ) )"),
                         [['complex', 'cplx', ['a', ['b', ['c', []]]]]])

        with self.assertRaises(ParseException):
            parse_kernel_string("cplx = a( b( c() ) )")

    def test_parse_examples(self):
        example1 = """
        length t0 = 6
        length d44 = 15
        length f0 = 15
        length d31 = 15
        length d25 = 15
        length Fluor25 = 15

        w44_31 = d44 t0 d31
        G31_25 = d31( t0( d25 + ) ) t0*
        R25    = d25( Fluor25 + ) t0*
        rep25  = d25 Fluor25
        #w31_f  = d31 t0 f0
        #w31_25 = d31 t0 d25
        #g31b   = t0* d31* t0*
        """
        output1 = [
            ['domain', 't0', '6'],
            ['domain', 'd44', '15'],
            ['domain', 'f0', '15'],
            ['domain', 'd31', '15'],
            ['domain', 'd25', '15'],
            ['domain', 'Fluor25', '15'],
            ['complex', 'w44_31', ['d44', 't0', 'd31']],
            ['complex', 'G31_25', ['d31', ['t0', ['d25', '+']], 't0*']],
            ['complex', 'R25', ['d25', ['Fluor25', '+'], 't0*']],
            ['complex', 'rep25', ['d25', 'Fluor25']]
        ]
        self.assertEqual(parse_kernel_string(example1), output1, 'seesaw example')

        self.assertEqual(parse_kernel_string("cplx = a( b( c( + ) ) d ) "),
                         [['complex', 'cplx', ['a', ['b', ['c', ['+']], 'd']]]])

        parse_kernel_string("cplx = a( b( c( + ) ) d ) @ constant 1e-7 M")
        parse_kernel_string("cplx = a( b( c( + ) ) d ) @ constant 1e-4 mM")
        parse_kernel_string("cplx = a( b( c( + ) ) d ) @ constant 0.1 uM")
        parse_kernel_string("cplx = a( b( c( + ) ) d ) @ initial 100 nM")
        parse_kernel_string("cplx = a( b( c( + ) ) d ) @ initial 1e5 pM")

        with self.assertRaises(ParseException):
            # whitespace between domains
            parse_kernel_string("cplx = a(b( c( + ) ) ) d ")

        with self.assertRaises(ParseException):
            # Closing domain must not be defined
            # NOTE: this behavior may change in the future?
            parse_kernel_string("cplx = a( b( c( + ) ) d) ")

        with self.assertRaises(ParseException):
            # Missing opening domain name
            parse_kernel_string("cplx = a( b( ( + ) ) ) ")

        with self.assertRaises(ParseException):
            # Unbalanced brackets
            parse_kernel_string("cplx = a( b( c( + ) ) d ")

        # NOTE: This test passes, ...
        # with self.assertRaises(ParseException):
        #  parse_kernel_string("cplx = a( b( c( + ) ) )d")


if __name__ == '__main__':
    unittest.main()