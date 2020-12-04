#
# tests/test_objectio.py
#   - copy and/or modify together with dsdobjects/objectio.py
#
import logging
logger = logging.getLogger('dsdobjects')
logger.setLevel(logging.INFO)
import unittest

import gc
from dsdobjects import SingletonError, clear_singletons
from dsdobjects.objectio import read_pil, read_pil_line, set_io_objects, clear_io_objects
from dsdobjects.base_classes import DomainS, StrandS, ComplexS, MacrostateS, ReactionS

SKIP = False

@unittest.skipIf(SKIP, "skipping tests.")
class TestReadLine(unittest.TestCase):
    def setUp(self):
        set_io_objects()

    def tearDown(self):
        clear_io_objects()

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
        assert repr(x) == 'ReactionS([condensed] A + I -> IA)'
        assert str(x) == '[condensed] A + I -> IA'
        assert x.reaction_string == 'reaction [condensed    = 1.66666e+06 /M/s ] A + I -> IA'

        self.assertEqual(list(x.products), [read_pil_line("macrostate IA = [IA]")])
        with self.assertRaises(SingletonError):
            read_pil_line("macrostate IA = [IA,A]")


@unittest.skipIf(SKIP, "skipping tests")
class TestReadFile(unittest.TestCase):
    def setUp(self):
        set_io_objects()

    def tearDown(self):
        clear_io_objects()

    def test_read_pil_01(self):
        out = read_pil(
        """
        # Domains (12) 
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
        assert len(out['con_reactions']) == 4
        A = out['macrostates']['A']
        I = out['macrostates']['I']
        IA = out['macrostates']['IA']
        assert ReactionS([A, I], [IA], 'condensed') in out['con_reactions']
        assert DomainS('a').sequence == 'NNNNNN'

    def test_complicated_input(self):
        out = read_pil(
        """
        # Domains (12) 
        sequence a = NNNNNN
        sequence b = NNNNNN 
        sequence c = NNNNNN
        sequence x = NNNNNN
        sequence y = NNNNNN 
        sequence z = NNNNNN 

        strand A = a x b y z* c* y* b* x*
        sup-sequence xby = x b y

        structure A = A : .(((..)))
        B = a xby( + ) b

        """)
        assert len(out['domains']) == 12
        assert len(out['strands']) == 2
        assert len(out['complexes']) == 2
        assert 'A' in out['strands']
        assert out['strands']['A'].structure is None
        assert 'xby' in out['strands']
        assert out['strands']['xby'].structure is None
        assert 'A' in out['complexes']
        assert out['complexes']['A'].kernel_string == 'a x( b( y( z* c* ) ) )'
        assert 'B' in out['complexes']
        assert list(out['complexes']['B'].sequence) == [DomainS('a'), DomainS('x'), DomainS('b'), DomainS('y'), '+', DomainS('y*'), DomainS('b*'), DomainS('x*'), DomainS('b')]
        assert out['complexes']['B'].kernel_string == 'a x( b( y( + ) ) ) b'

    def test_macrostate_name_bug(self):
        out = read_pil("""
            # Domains (8) 
            length d2 = 15
            length d3 = 15
            length d4 = 15
            length t1 = 5

            # Resting complexes (32) 
            A = t1 d2
            f1 = d2( t1( + d4( + d3( t1( + d3( + ) ) ) ) ) ) t1* @constant 100 nM
            f2 = t1 d4 d3 @constant 100 nM
            f3 = d2( + t1( d3( + t1* ) ) ) @constant 100 nM

            # ...
            e36 = t1 d3
            e499 = t1( d4 d3( + t1* ) ) d2*( + t1 )
            e501 = t1 d4 d3( + t1* ) t1* d2*( + t1 )
            e1027 = t1* d3*( t1*( d2*( + t1 ) + ) )
            # ...
            macrostate f2 = [f2]
            macrostate e36 = [e36]
            macrostate e499 = [e501, e499]
            macrostate e1027 = [e1027]
            # ...
            reaction [condensed    = 6.57478e-05 /nM/s ] e499 + e36 -> e1027 + f2
            # ...
            """)
        pass

if __name__ == '__main__':
    unittest.main()

