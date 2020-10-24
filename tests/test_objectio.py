#
# tests/test_objectio.py
#   - copy and/or modify together with dsdobjects/objectio.py
#
import logging
logger = logging.getLogger('dsdobjects')
logger.setLevel(logging.INFO)
import unittest

from dsdobjects import SingletonError, clear_singletons
from dsdobjects.objectio import read_pil, read_pil_line, set_prototypes
from dsdobjects.base_classes import DomainS, ComplexS, MacrostateS, ReactionS

SKIP = False

@unittest.skipIf(SKIP, "skipping tests.")
class Test_IO(unittest.TestCase):
    def setUp(self):
        set_prototypes()

    def tearDown(self):
        clear_singletons(DomainS)
        clear_singletons(ComplexS)
        clear_singletons(MacrostateS)
        clear_singletons(ReactionS)

    def test_read_pil_line_01(self):
        x = read_pil_line("length d5 = 15")

        assert isinstance(x, DomainS)
        assert x.name == 'd5'
        assert x.length == 15

        with self.assertRaises(SingletonError):
            x = read_pil_line("length d5 = 12")

        x = read_pil_line("sequence t1 = NCGGA")
        assert isinstance(x, DomainS)
        assert x.name == 't1'
        assert x.length == 5
        assert x.sequence == 'NCGGA'
        assert (~x).sequence == 'TCCGN'

    def test_read_pil_line_02(self):
        doms = []
        doms.append(read_pil_line("length a = 6"))
        doms.append(read_pil_line("length b = 6"))
        doms.append(read_pil_line("length c = 6"))
        doms.append(read_pil_line("length x = 6"))
        doms.append(read_pil_line("length y = 6"))
        doms.append(read_pil_line("length z = 6"))
        for d in list(doms):
            doms.append(~d)

        A = read_pil_line("A = a x( b( y( z* c* ) ) )")
        I = read_pil_line("I = y* b* x* a*")
        IA = read_pil_line("IA = a( x( b( y( z* c* y* b* x* + ) ) ) )")
        assert isinstance(IA, ComplexS)
        assert IA.name == 'IA'

    def test_read_pil_line_03(self):
        doms = []
        doms.append(read_pil_line("sequence a = NNNNNN"))
        doms.append(read_pil_line("sequence b = NNNNNN"))
        doms.append(read_pil_line("sequence c = NNNNNN"))
        doms.append(read_pil_line("sequence x = NNNNNN"))
        doms.append(read_pil_line("sequence y = NNNNNN"))
        doms.append(read_pil_line("sequence z = NNNNNN"))
        for d in list(doms):
            doms.append(~d)
        A = read_pil_line("A = a x( b( y( z* c* ) ) )")
        I = read_pil_line("I = y* b* x* a*")
        IA = read_pil_line("IA = a( x( b( y( z* c* y* b* x* + ) ) ) )")
        assert isinstance(IA, ComplexS)
        assert IA.name == 'IA'

        A = read_pil_line("macrostate A = [A]")
        I = read_pil_line("macrostate I = [I]")
        IA = read_pil_line("macrostate IA = [IA]")

        assert isinstance(IA, MacrostateS)
        assert IA.name == 'IA'
        assert list(IA.complexes) == [read_pil_line("IA = a( x( b( y( z* c* y* b* x* + ) ) ) )")]

        x = read_pil_line("reaction [condensed      =  1.66666e+06 /M/s ] A + I -> IA")
        assert isinstance(x, ReactionS)
        assert repr(x) == "ReactionS(reaction [condensed = 1666660 /M/s] A + I -> IA)"
        assert str(x) == "reaction [condensed = 1666660 /M/s] A + I -> IA"

        self.assertEqual(x.products, [read_pil_line("macrostate IA = [IA]")])
        with self.assertRaises(SingletonError):
            read_pil_line("macrostate IA = [IA,A]")


@unittest.skipIf(SKIP, "skipping tests")
class Test_Reaction(unittest.TestCase):
    def setUp(self):
        set_prototypes()

    def tearDown(self):
        clear_singletons(DomainS)
        clear_singletons(ComplexS)
        clear_singletons(MacrostateS)
        clear_singletons(ReactionS)

    def test_ptreact(self):
        dom, clx, rms, det, con = read_pil(
        """
        # Domains (12) 
        length a = 6
        length b = 6
        length c = 6
        length x = 6
        length y = 6
        length z = 6

        sequence a = NNNNNN
        sequence b = NNNNNN 
        sequence c = NNNNNN
        sequence x = NNNNNN
        sequence y = NNNNNN 
        sequence z = NNNNNN 

        # Resting complexes (8) 
        A = a x( b( y( z* c* ) ) ) @i 1e-08 M
        B = b y( c( z( x* a* ) ) ) @i 1e-08 M
        C = c z( a( x( y* b* ) ) ) @i 1e-08 M
        I = y* b* x* a* @i 1e-08 M

        IA = a( x( b( y( z* c* y* b* x* + ) ) ) ) @i 0.0 M
        IAB = y*( b*( x*( a*( + ) ) ) ) z*( c*( y*( b*( x* + ) ) ) ) x* a* z* c* y* @i 0.0 M
        IABC = y*( b*( x*( a*( + ) ) ) ) z*( c*( y*( b*( x* + ) ) ) ) x*( a*( z*( c*( y* + ) ) ) ) y* b* x* a* z* @i 0.0 M
        ABC = a( x( b( y( z*( c*( y*( b*( x* + ) ) ) ) x*( a*( z*( c*( y* + ) ) ) ) ) ) ) ) z* @i 0.0 M

        # Resting macrostates (8) 
        macrostate A = [A]
        macrostate B = [B]
        macrostate C = [C]
        macrostate I = [I]
        macrostate IA = [IA]
        macrostate IAB = [IAB]
        macrostate IABC = [IABC]
        macrostate ABC = [ABC]

        # Condensed reactions (10) 
        reaction [condensed      =  1.66666e+06 /M/s ] A + I -> IA
        reaction [condensed      =  1.66666e+06 /M/s ] IA + B -> IAB
        reaction [condensed      =  1.66646e+06 /M/s ] IAB + C -> IABC
        reaction [condensed      =    0.0261637 /s   ] IABC -> ABC + I
        """)

        # A preliminary interface to start testing prototype functions.
        rxn = con[0]
        a, b, c = dom['a'], dom['b'], dom['c']
        x, y, z = dom['x'], dom['y'], dom['z']
        A, B, I = rms['A'], rms['B'], rms['I']

        assert rxn.reactants == [A, I]
        #so, pt1, pt2 = rxn.ptreact()
        #assert so.sequence == [a, x, b, y, ~z, ~c, ~y, ~b, ~x, '+', ~y, ~b, ~x, ~a]
        #assert pt1 == [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)], [None, None, None, None]]
        #assert pt2 == [[(1, 3), (1, 2), (1, 1), (1, 0), None, None, None, None, None], [(0, 3), (0, 2), (0, 1), (0, 0)]]

if __name__ == '__main__':
    unittest.main()

