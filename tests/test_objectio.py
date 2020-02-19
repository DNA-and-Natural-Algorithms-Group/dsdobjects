import unittest

from dsdobjects import LogicDomain, Domain, Complex, Reaction, Macrostate, StrandOrder
from dsdobjects import DSDObjectsError
from dsdobjects.core import clear_memory
from dsdobjects.objectio import read_pil, read_pil_line, set_prototypes


class Test_IO(unittest.TestCase):
    def setUp(self):
        set_prototypes()

    def tearDown(self):
        clear_memory()

    def test_read_pil_line(self):
        x = read_pil_line("length d5 = 15")

        assert isinstance(x, LogicDomain)
        assert x.name == 'd5'
        assert x.length == 15
        assert x.dtype == 'long'

        with self.assertRaises(DSDObjectsError):
            x = read_pil_line("length d5 = 12")

        x = read_pil_line("sequence t1 = NCGGA")

        assert isinstance(x, Domain)
        assert x.name == 't1'
        assert x.length == 5
        assert x.dtype == LogicDomain('t1', length=5)

        x = read_pil_line("length a = 6")  
        x = read_pil_line("length b = 6")
        x = read_pil_line("length c = 6")
        x = read_pil_line("length x = 6")
        x = read_pil_line("length y = 6")
        x = read_pil_line("length z = 6")

        x = read_pil_line("sequence a = NNNNNN")  
        x = read_pil_line("sequence b = NNNNNN")
        x = read_pil_line("sequence c = NNNNNN")
        x = read_pil_line("sequence x = NNNNNN")
        x = read_pil_line("sequence y = NNNNNN")
        x = read_pil_line("sequence z = NNNNNN")

        x = read_pil_line("A = a x( b( y( z* c* ) ) )")
        x = read_pil_line("I = y* b* x* a*")
        x = read_pil_line("IA = a( x( b( y( z* c* y* b* x* + ) ) ) )")
        assert isinstance(x, Complex)
        assert x.name == 'IA'
        assert x.size == 2

        x = read_pil_line("macrostate A = [A]")
        x = read_pil_line("macrostate I = [I]")
        x = read_pil_line("macrostate IA = [IA]")
        assert isinstance(x, Macrostate)
        assert x.name == 'IA'
        assert x.complexes == [read_pil_line("IA = a( x( b( y( z* c* y* b* x* + ) ) ) )")]

        x = read_pil_line("reaction [condensed      =  1.66666e+06 /M/s ] A + I -> IA")
        assert isinstance(x, Reaction)
        self.assertEqual(str(x), 'A + I -> IA')
        self.assertEqual(x.full_string(), "[condensed      =  1.66666e+06 /M/s ] A + I -> IA")
        self.assertEqual(x.products, [read_pil_line("macrostate IA = [IA]")])

        with self.assertRaises(DSDObjectsError):
            read_pil_line("macrostate IA = [IA,A]")

if __name__ == '__main__':
    unittest.main()

