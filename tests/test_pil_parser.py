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

        with self.assertRaises(ParseException):
            parse_pil_string(" state e4 = e4")

    def test_not_implemented(self):
        # there is currently no support of additional arguments
        # for complexes, e.g. this: [1nt]
        with self.assertRaises(ParseException):
            parse_pil_string(" structure [1nt] AB = A B : .(((+))) ")
        with self.assertRaises(ParseException):
            parse_pil_string(" structure [1nt] AB = A + B : ......((((((((((((((((((((((((((((((((((+)))))))))))))))))))))))))))))))))) ")

    def test_valid_inputs(self):
        # use print statements to see return values
        out = parse_pil_string(" domain a : 6                    ")
        out = parse_pil_string(" domain a : short                ")
        self.assertEqual(out, [['dl-domain', 'a', 'short']])
        out = parse_pil_string(" sequence a : 18                 ")
        self.assertEqual(out, [['dl-domain', 'a', '18']])
        out = parse_pil_string(" sequence a1 = CTAGA : 6         ")
        self.assertEqual(out, [['sl-domain', 'a1', 'CTAGA', '6']])
        out = parse_pil_string(" sequence a1 = CTAGA             ")
        self.assertEqual(out, [['sl-domain', 'a1', 'CTAGA']])
        out = parse_pil_string(" sequence 1 = CTA : 6            ")
        self.assertEqual(out, [['sl-domain', '1', 'CTA', '6']])

        # composite-domain
        out = parse_pil_string(" sup-sequence q = a b-seq z : 20 ")
        self.assertEqual(out, [['composite-domain', 'q', ['a', 'b-seq', 'z'], '20']])
        out = parse_pil_string(" strand q = a b-seq z : 20       ")
        self.assertEqual(out, [['composite-domain', 'q', ['a', 'b-seq', 'z'], '20']])
        out = parse_pil_string(" strand q = a b-seq z            ")
        self.assertEqual(out, [['composite-domain', 'q', ['a', 'b-seq', 'z']]])
        out = parse_pil_string(" strand I : y* b* x* a*          ")
        self.assertEqual(out, [['composite-domain', 'I', ['y*', 'b*', 'x*', 'a*']]])
        out = parse_pil_string(" structure AB = A B : .(((+)))   ")
        self.assertEqual(out, [['strand-complex', 'AB', ['A', 'B'], '.(((+)))   ']])
        out = parse_pil_string(" structure AB = A + B : .(((+))) ")
        self.assertEqual(out, [['strand-complex', 'AB', ['A', 'B'], '.(((+))) ']])
        out = parse_pil_string(" kinetic 4 + C1 -> 7             ")
        self.assertEqual(out, [['reaction', [], ['4', 'C1'], ['7']]])
        out = parse_pil_string(" kinetic [   876687.69 /M/s] I3 + SP -> C2 + Cat ")
        self.assertEqual(out, [['reaction', [[], ['876687.69']], ['I3', 'SP'], ['C2', 'Cat']]])
        out = parse_pil_string(" kinetic [   1667015.4 /M/s] I3 + C1 -> W + Cat + OP ")
        self.assertEqual(out, [['reaction', [[], ['1667015.4']], ['I3', 'C1'], ['W', 'Cat', 'OP']]])
        out = parse_pil_string(" state e4 = [e4]")
        self.assertEqual(out, [['resting-set', 'e4', ['e4']]])
        out = parse_pil_string(" state e4 = [e4, e5]")
        self.assertEqual(out, [['resting-set', 'e4', ['e4', 'e5']]])

        # reaction
        out = parse_pil_string(" reaction [branch-3way =  0.733333 /s   ] e71 -> e11 ")
        self.assertEqual(out, [['reaction', [['branch-3way'], ['0.733333']], ['e71'], ['e11']]])
        out = parse_pil_string(" reaction [bind21      =   4.5e+06 /M/s ] e4 + G1bot -> e13")
        self.assertEqual(out, [['reaction', [['bind21'], ['4.5e+06']], ['e4', 'G1bot'], ['e13']]])

        # strand-complex
        out = parse_pil_string(" complex IABC : I A B C (((( + ))))((((. + ))))((((. + ))))..... ")
        self.assertEqual(out, [['strand-complex', 'IABC', ['I', 'A', 'B', 'C'], '(((( + ))))((((. + ))))((((. + ))))..... ']])
        out = parse_pil_string(" complex I : I ....              ") 
        self.assertEqual(out, [['strand-complex', 'I', ['I'], '....              ']])

        # kernel-complex
        out = parse_pil_string(" e10 = 2( 3 + 3( 4( + ) ) ) 1*( + ) 2  @ initial 0 nM")
        self.assertEqual(out, [['kernel-complex', 'e10', ['2', ['3', '+', '3', ['4', ['+']]], '1*', ['+'], '2'], ['initial', '0', 'nM']]])
        out = parse_pil_string(" C = 1 2 3( ) + 4 ")
        self.assertEqual(out, [['kernel-complex', 'C', ['1', '2', '3', [], '+', '4']]])
        out = parse_pil_string(" C = 1 2 3() + 4 ")
        self.assertEqual(out, [['kernel-complex', 'C', ['1', '2', '3', [], '+', '4']]])
        out = parse_pil_string(" C = 1 2(3(+))")
        self.assertEqual(out, [['kernel-complex', 'C', ['1', '2', ['3', ['+']]]]])

        out = parse_pil_string(" fuel2 = 2b( 3a( 3b( 3c( 3d( 4a( 4b 4c + ) ) ) ) ) ) 2a*         ")
        self.assertEqual(out, [['kernel-complex', 'fuel2', ['2b', ['3a', ['3b', ['3c', ['3d', ['4a', ['4b', '4c', '+']]]]]], '2a*']]])
        out = parse_pil_string(" fuel2 = 2b( 3a( 3b( 3c( 3d( 4a( 4b 4c + )))))) 2a*              ")
        self.assertEqual(out, [['kernel-complex', 'fuel2', ['2b', ['3a', ['3b', ['3c', ['3d', ['4a', ['4b', '4c', '+']]]]]], '2a*']]])
        out = parse_pil_string(" fuel2 = 2b(3a( 3b( 3c( 3d( 4a( 4b 4c + )))))) 2a*               ")
        self.assertEqual(out, [['kernel-complex', 'fuel2', ['2b', ['3a', ['3b', ['3c', ['3d', ['4a', ['4b', '4c', '+']]]]]], '2a*']]])
        out = parse_pil_string(" fuel2 = 2b(3a( 3b( 3c( 3d( 4a( )))))) 2a*                       ")
        self.assertEqual(out, [['kernel-complex', 'fuel2', ['2b', ['3a', ['3b', ['3c', ['3d', ['4a', []]]]]], '2a*']]])

        out = parse_pil_string("""
            complex IABC : 
            I A B C 
            (((( + ))))((((. + ))))((((. + ))))..... 
            """)
        self.assertEqual(out, [['strand-complex', 'IABC', ['I', 'A', 'B', 'C'], '(((( + ))))((((. + ))))((((. + ))))..... ']])

if __name__ == '__main__':
    unittest.main()
