import unittest
from pyparsing import ParseException
from dsdobjects.parser import parse_kernel_file, parse_kernel_string
from dsdobjects.parser import parse_pil_file, parse_pil_string

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class TestPILparser(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_broken_inputs(self):
        # from peppercornenumerator/tests
        with self.assertRaises(ParseException):
            parse_pil_string(""" sequence a = CT1 : 6""")

        with self.assertRaises(ParseException):
            parse_pil_string(""" strand = CT1 : 6""")

    def test_not_implemented(self):
        # there is currently no support of werid [1nt] keyword:
        with self.assertRaises(ParseException):
            parse_pil_string(" structure [1nt] AB = A B : .(((+))) ")
        with self.assertRaises(ParseException):
            parse_pil_string(" structure [1nt] AB = A + B : ......((((((((((((((((((((((((((((((((((+)))))))))))))))))))))))))))))))))) ")

    def test_valid_inputs(self):
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
