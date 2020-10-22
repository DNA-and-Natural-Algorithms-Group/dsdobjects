#
# tests for dsdobjects.dsdparser.seesaw_parser.py
#
import unittest
from pyparsing import ParseException
from dsdobjects.dsdparser import parse_seesaw_string, parse_seesaw_file

class TestSeeSawParser(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_io(self):
        out = parse_seesaw_string("INPUT(1) = w[44, 31]")
        self.assertEqual(out, [['INPUT', ['1'], ['w', ['44', '31']]]])

        out = parse_seesaw_string("OUTPUT(1) = w[44, 31]")
        self.assertEqual(out, [['OUTPUT', ['1'], ['w', ['44', '31']]]])

        out = parse_seesaw_string("OUTPUT(2) = Fluor[25]")
        self.assertEqual(out, [['OUTPUT', ['2'], ['Fluor', '25']]])

        out = parse_seesaw_string("INPUT(x1) = w[44, 31]")
        self.assertEqual(out, [['INPUT', ['x1'], ['w', ['44', '31']]]])

        with self.assertRaises(ParseException):
            print(parse_seesaw_string("INPUT(1) = Fluor[25]"))

        with self.assertRaises(ParseException):
            print(parse_seesaw_string("OUTPUT(21x) = Fluor[25]"))

    def test_seesaw(self):
        out = parse_seesaw_string("seesaw[31, {44}, {25, f}]")
        self.assertEqual(out, [['seesaw', ['31', ['44'], ['25', 'f']]]])

        with self.assertRaises(ParseException):
            print(parse_seesaw_string("seesaw[31, 19, {44}, {25, f}]"))

        with self.assertRaises(ParseException):
            print(parse_seesaw_string("seesaw[31, 44, {25, f}]"))


    def test_gates(self):
        out = parse_seesaw_string("seesawAND[1, 44, {3, 9, 17}, {13, 19, 117}]")
        self.assertEqual(out, [['seesawAND', ['1', '44', ['3', '9', '17'], ['13', '19', '117']]]])

        out = parse_seesaw_string("seesawOR[1, 44, {3, 9, 17}, {13, 19, 117}]")
        self.assertEqual(out, [['seesawOR', ['1', '44', ['3', '9', '17'], ['13', '19', '117']]]])

        out = parse_seesaw_string("inputfanout[1, 44, {3, 9,17}]")
        self.assertEqual(out, [['inputfanout', ['1', '44', ['3', '9', '17']]]])

        out = parse_seesaw_string("reporter[13, 44]")
        self.assertEqual(out, [['reporter', ['13', '44']]])

        with self.assertRaises(ParseException):
            print(parse_seesaw_string("seesawOR[1, {3, 9, 17}, {13, 19, 117}]"))
        
        with self.assertRaises(ParseException):
            print(parse_seesaw_string("inputfanout[1, {44}, {3, 9,17}]"))

    def test_concentrations(self):
        #wires
        out = parse_seesaw_string("conc[w[31, 44], 1*c]")
        self.assertEqual(out, [['conc', ['w', ['31', '44']], '1']])
        out = parse_seesaw_string("conc[w[31, f], 2*c]")
        self.assertEqual(out, [['conc', ['w', ['31', 'f']], '2']])

        #gates
        out = parse_seesaw_string("conc[g[w[31, 44], 13], 1*c]")
        self.assertEqual(out, [['conc', ['g', [['w', ['31', '44']], '13']], '1']])
        out = parse_seesaw_string("conc[g[13, w[31, 44]], 1*c]")
        self.assertEqual(out, [['conc', ['g', ['13', ['w', ['31', '44']]]], '1']])

        #thresholds
        out = parse_seesaw_string("conc[th[w[31, 44], 13], 1*c]")
        self.assertEqual(out, [['conc', ['th', [['w', ['31', '44']], '13']], '1']])
        out = parse_seesaw_string("conc[th[13, w[31, 44]], 1*c]")
        self.assertEqual(out, [['conc', ['th', ['13', ['w', ['31', '44']]]], '1']])

        #non-integer concentration
        out = parse_seesaw_string("conc[th[13, w[31, 44]], 0.5*c]")
        self.assertEqual(out, [['conc', ['th', ['13', ['w', ['31', '44']]]], '0.5']])

        #negative concentration
        with self.assertRaises(ParseException):
            print(parse_seesaw_string("conc[th[13, w[31, 44]], -0.5*c]"))

        with self.assertRaises(ParseException):
            print(parse_seesaw_string("conc[th[13, w[31, 44]], -1*c]"))

    def test_parse_examples(self):
        example = """
        # Example catalytic circuit (Lab 6)

        INPUT(x) = w[44, 31]    # x
        OUTPUT(y) = Fluor[25]   # y

        seesaw[31, {44}, {25, f}]
        conc[w[31,f],2*c]
        conc[g[31, w[31,25]],1*c]
        reporter[25, 31]
        """
        output = [['INPUT', ['x'],  ['w', ['44', '31']]], 
                  ['OUTPUT', ['y'], ['Fluor', '25']], 
                  ['seesaw', ['31', ['44'], ['25', 'f']]], 
                  ['conc', ['w', ['31', 'f']], '2'],
                  ['conc', ['g', ['31', ['w', ['31', '25']]]], '1'], 
                  ['reporter', ['25', '31']]]
        self.assertEqual(parse_seesaw_string(example), output)

        example = """
        # a seesaw analog circuit that generates a pulse

        INPUT(1) = w[1,2]    
        OUTPUT(2) = w[2,3]  

        seesaw[2,{1,4},{3}]
        seesaw[4,{6},{5,2}]

        conc[g[w[4,2],2],2*c]
        conc[g[2,w[2,3]],1*c]
        conc[g[w[6,4],4],1*c]
        conc[g[4,w[4,5]],100*c]
        """
        output = [['INPUT', ['1'], ['w', ['1', '2']]], 
                  ['OUTPUT', ['2'], ['w', ['2', '3']]], 

                  ['seesaw', ['2', ['1', '4'], ['3']]], 
                  ['seesaw', ['4', ['6'], ['5', '2']]], 

                  ['conc', ['g', [['w', ['4', '2']], '2']], '2'], 
                  ['conc', ['g', ['2', ['w', ['2', '3']]]], '1'], 
                  ['conc', ['g', [['w', ['6', '4']], '4']], '1'], 
                  ['conc', ['g', ['4', ['w', ['4', '5']]]], '100']]
        self.assertEqual(parse_seesaw_string(example), output)


if __name__ == '__main__':
  unittest.main()

