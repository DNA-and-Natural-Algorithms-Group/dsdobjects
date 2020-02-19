import unittest

from dsdobjects import LogicDomain, Domain, Complex, Reaction, Macrostate, StrandOrder
from dsdobjects import DSDObjectsError
from dsdobjects.core import clear_memory
from dsdobjects.objectio import read_pil, set_prototypes

set_prototypes()
SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class Test_Reaction(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

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
        so,pt1,pt2 = rxn.ptreact()

        assert so.sequence == [a, x, b, y, ~z, ~c, ~y, ~b, ~x, '+', ~y, ~b, ~x, ~a]
        assert pt1 == [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)], [None, None, None, None]]
        assert pt2 == [[(1, 3), (1, 2), (1, 1), (1, 0), None, None, None, None, None], [(0, 3), (0, 2), (0, 1), (0, 0)]]


if __name__ == '__main__':
    unittest.main()

