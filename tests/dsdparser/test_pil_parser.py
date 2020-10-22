#
# tests for dsdobjects.dsdparser.pil_parser.py
#
import unittest
from pyparsing import ParseException

from dsdobjects.dsdparser import parse_pil_file, parse_pil_string

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

    def test_dl_domains(self):
        out = parse_pil_string(" length a : short                ")
        self.assertEqual(out, [['dl-domain', 'a', 'short']])
        out = parse_pil_string(" domain a : 6                    ")
        self.assertEqual(out, [['dl-domain', 'a', '6']])
        out = parse_pil_string(" domain a : short                ")
        self.assertEqual(out, [['dl-domain', 'a', 'short']])
        out = parse_pil_string(" sequence a : 18                 ")
        self.assertEqual(out, [['dl-domain', 'a', '18']])

    def test_sl_domains(self):
        out = parse_pil_string(" sequence a1 = CTAGA : 6         ")
        self.assertEqual(out, [['sl-domain', 'a1', 'CTAGA', '6']])
        out = parse_pil_string(" sequence a1 = CTAGA             ")
        self.assertEqual(out, [['sl-domain', 'a1', 'CTAGA']])
        out = parse_pil_string(" sequence 1 = CTA : 6            ")
        self.assertEqual(out, [['sl-domain', '1', 'CTA', '6']])

    def test_composite_domains(self):
        out = parse_pil_string(" sup-sequence q = a b-seq z : 20 ")
        self.assertEqual(out, [['composite-domain', 'q', ['a', 'b-seq', 'z'], '20']])
        out = parse_pil_string(" strand q = a b-seq z : 20       ")
        self.assertEqual(out, [['composite-domain', 'q', ['a', 'b-seq', 'z'], '20']])
        out = parse_pil_string(" strand q = a b-seq z            ")
        self.assertEqual(out, [['composite-domain', 'q', ['a', 'b-seq', 'z']]])
        out = parse_pil_string(" strand I : y* b* x* a*          ")
        self.assertEqual(out, [['composite-domain', 'I', ['y*', 'b*', 'x*', 'a*']]])

    def test_strand_complex(self):
        out = parse_pil_string(" structure AB = A B : .(((+)))   ")
        self.assertEqual(out, [['strand-complex', 'AB', ['A', 'B'], '.(((+)))   ']])

        out = parse_pil_string(" structure AB = A + B : .(((+))) ")
        self.assertEqual(out, [['strand-complex', 'AB', ['A', 'B'], '.(((+))) ']])

        out = parse_pil_string(" complex IABC : I A B C (((( + ))))((((. + ))))((((. + ))))..... ")
        self.assertEqual(out, [['strand-complex', 'IABC', ['I', 'A', 'B', 'C'], '(((( + ))))((((. + ))))((((. + ))))..... ']])

        out = parse_pil_string(" complex I : I ....              ") 
        self.assertEqual(out, [['strand-complex', 'I', ['I'], '....              ']])

        out = parse_pil_string("""
            complex IABC : 
            I A B C 
            (((( + ))))((((. + ))))((((. + ))))..... 
            """)
        self.assertEqual(out, [['strand-complex', 'IABC', ['I', 'A', 'B', 'C'], '(((( + ))))((((. + ))))((((. + ))))..... ']])

    def test_reactions(self):
        out = parse_pil_string(" kinetic 4 + C1 -> 7             ")
        self.assertEqual(out, [['reaction', [], ['4', 'C1'], ['7']]])

        out = parse_pil_string(" kinetic [   876687.69 /M/s] I3 + SP -> C2 + Cat ")
        self.assertEqual(out, [['reaction', [[], ['876687.69'], ['/M/s']], ['I3', 'SP'], ['C2', 'Cat']]])

        out = parse_pil_string(" kinetic [   1667015.4 /M/s] I3 + C1 -> W + Cat + OP ")
        self.assertEqual(out, [['reaction', [[], ['1667015.4'], ['/M/s']], ['I3', 'C1'], ['W', 'Cat', 'OP']]])

        out = parse_pil_string(" reaction [branch-3way =  0.733333 /s   ] e71 -> e11 ")
        self.assertEqual(out, [['reaction', [['branch-3way'], ['0.733333'], ['/s']], ['e71'], ['e11']]])

        out = parse_pil_string(" reaction [bind21      =   4.5e+06 /M/s ] e4 + G1bot -> e13")
        self.assertEqual(out, [['reaction', [['bind21'], ['4.5e+06'], ['/M/s']], ['e4', 'G1bot'], ['e13']]])

        out = parse_pil_string(" reaction [bind21      =   4.5e+06 /nM/h ] e4 + G1bot -> e13")
        self.assertEqual(out, [['reaction', [['bind21'], ['4.5e+06'], ['/nM/h']], ['e4', 'G1bot'], ['e13']]])

        out = parse_pil_string("reaction [k1 =  1.41076e+07 +/-  1.47099e+06 /M/s] A + B -> A_B")
        self.assertEqual(out, [['reaction', [['k1'], ['1.41076e+07', '1.47099e+06'], ['/M/s']], ['A', 'B'], ['A_B']]])

        out = parse_pil_string("reaction [k1 =  1.41076e+07 +/-  inf /M/s] A + B -> A_B")
        self.assertEqual(out, [['reaction', [['k1'], ['1.41076e+07', 'inf'], ['/M/s']], ['A', 'B'], ['A_B']]])

    def test_resting_macrostate(self):
        out = parse_pil_string(" state e4 = [e4]")
        self.assertEqual(out, [['resting-macrostate', 'e4', ['e4']]])
        out = parse_pil_string(" state e4 = [e4, e5]")
        self.assertEqual(out, [['resting-macrostate', 'e4', ['e4', 'e5']]])

    def test_kernel_complex(self):
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

    def test_complex_concentrations(self):
        out = parse_pil_string("cplx = a( b( c( + ) ) d ) @ constant 1e-7 M")
        self.assertEqual(out, [['kernel-complex', 'cplx', ['a', ['b', ['c', ['+']], 'd']], ['constant', '1e-7', 'M']]])
        out = parse_pil_string("cplx = a( b( c( + ) ) d ) @ constant 1e-4 mM")
        self.assertEqual(out, [['kernel-complex', 'cplx', ['a', ['b', ['c', ['+']], 'd']], ['constant', '1e-4', 'mM']]])
        out = parse_pil_string("cplx = a( b( c( + ) ) d ) @ constant 0.1 uM")
        self.assertEqual(out, [['kernel-complex', 'cplx', ['a', ['b', ['c', ['+']], 'd']], ['constant', '0.1', 'uM']]])
        out = parse_pil_string("cplx = a( b( c( + ) ) d ) @ initial 100 nM")
        self.assertEqual(out, [['kernel-complex', 'cplx', ['a', ['b', ['c', ['+']], 'd']], ['initial', '100', 'nM']]])
        out = parse_pil_string("cplx = a( b( c( + ) ) d ) @ initial 1e5 pM")
        self.assertEqual(out, [['kernel-complex', 'cplx', ['a', ['b', ['c', ['+']], 'd']], ['initial', '1e5', 'pM']]])

if __name__ == '__main__':
    unittest.main()
