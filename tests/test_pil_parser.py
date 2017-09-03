# -*- coding: utf-8 -*-
#
# tests for dsdobjects.parser.pil_extended_format
#
# Written by Stefan Badelt (badelt@caltech.edu)
#

import unittest
from pyparsing import ParseException

from dsdobjects.parser import parse_pil_file, parse_pil_string

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class TestPILparser(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_broken_inputs(self):
        # A digit in the sequence ...
        with self.assertRaises(ParseException):
            parse_pil_string(""" sequence a = CT1 : 6""")

        # A strand with no name ...
        with self.assertRaises(ParseException):
            parse_pil_string(""" strand = CT1 : 6""")

    def test_not_implemented(self):
        # there is currently no support of additional arguments
        # for complexes, e.g. this: [1nt]
        with self.assertRaises(ParseException):
            parse_pil_string(" structure [1nt] AB = A B : .(((+))) ")
        with self.assertRaises(ParseException):
            parse_pil_string(" structure [1nt] AB = A + B : ......((((((((((((((((((((((((((((((((((+)))))))))))))))))))))))))))))))))) ")

    def test_valid_inputs(self):
        # use print statements to see return values
        parse_pil_string(" domain a : 6                    ")
        parse_pil_string(" domain a : short                ")
        parse_pil_string(" sequence a : 18                 ")
        parse_pil_string(" sequence a1 = CTAGA : 6         ")
        parse_pil_string(" sequence a1 = CTAGA             ")
        parse_pil_string(" sequence 1 = CTA : 6            ")
        parse_pil_string(" sup-sequence q = a b-seq z : 20 ")
        parse_pil_string(" strand q = a b-seq z : 20       ")
        parse_pil_string(" strand q = a b-seq z            ")
        parse_pil_string(" strand I : y* b* x* a*          ")
        parse_pil_string(" structure AB = A B : .(((+)))   ")
        parse_pil_string(" structure AB = A + B : .(((+))) ")
        parse_pil_string(" complex I : I ....              ")
        parse_pil_string(" kinetic 4 + C1 -> 7             ")
        parse_pil_string(" kinetic [   876687.69 /M/s] I3 + SP -> C2 + Cat ")
        parse_pil_string(" kinetic [   1667015.4 /M/s] I3 + C1 -> W + Cat + OP ")
        parse_pil_string(" reaction [branch-3way =  0.733333 /s   ] e71 -> e11 ")
        parse_pil_string(" reaction [bind21      =   4.5e+06 /M/s ] e4 + G1bot -> e13")
        parse_pil_string(" complex IABC : I A B C (((( + ))))((((. + ))))((((. + ))))..... ")
        parse_pil_string(" fuel2 = 2b( 3a( 3b( 3c( 3d( 4a( 4b 4c + ) ) ) ) ) ) 2a*         ")
        parse_pil_string(" fuel2 = 2b( 3a( 3b( 3c( 3d( 4a( 4b 4c + )))))) 2a*              ")
        parse_pil_string(" fuel2 = 2b(3a( 3b( 3c( 3d( 4a( 4b 4c + )))))) 2a*               ")
        parse_pil_string(" fuel2 = 2b(3a( 3b( 3c( 3d( 4a( )))))) 2a*                       ")
        parse_pil_string("""
            complex IABC : 
            I A B C 
            (((( + ))))((((. + ))))((((. + ))))..... 
            """)

        # Outputs:
        # [['dl-domain', 'a', '6']]
        # [['sl-domain', 'a1', 'CTAGA', '6']]
        # [['sl-domain', 'a1', 'CTAGA']]
        # [['sl-domain', '1', 'CTA', '6']]

        # [['composite-domain', 'q', ['a', 'b-seq', 'z'], '20']]
        # [['composite-domain', 'q', ['a', 'b-seq', 'z'], '20']]
        # [['composite-domain', 'q', ['a', 'b-seq', 'z']]]
        # [['composite-domain', 'I', ['y*', 'b*', 'x*', 'a*']]]

        # [['strand-complex', 'AB', ['A', 'B'], '.(((+)))   ']]
        # [['strand-complex', 'AB', ['A', 'B'], '.(((+))) ']]
        # [['strand-complex', 'I', ['I'], '....              ']]
        # [['strand-complex', 'IABC', ['I', 'A', 'B', 'C'], '(((( + ))))((((. + ))))((((. + ))))..... ']]
        # [['strand-complex', 'IABC', ['I', 'A', 'B', 'C'], '(((( + ))))((((. + ))))((((. + ))))..... ']]

if __name__ == '__main__':
    unittest.main()
