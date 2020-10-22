#
# Parser module for Seesaw compiler specification of DSD cirucits.
# http://www.qianlab.caltech.edu/SeesawCompiler/
#
# dsdobjects.dsdparser.seesaw_parser.py
#   - copy and/or modify together with tests/dsdparser/test_seesaw_parser.py
#
from pyparsing import (Word, Literal, Group, Suppress, Optional, ZeroOrMore, Combine, 
    OneOrMore, alphas, alphanums, nums, delimitedList, StringStart, StringEnd, 
    LineEnd, pythonStyleComment, ParseElementEnhance)

def ssw_document_setup():
  crn_DWC = "".join(
      [x for x in ParseElementEnhance.DEFAULT_WHITE_CHARS if x != "\n"])
  ParseElementEnhance.setDefaultWhitespaceChars(crn_DWC)
 
  def T(x, tag):
    def TPA(tag):
      return lambda s, l, t: [tag] + t.asList()
    return x.setParseAction(TPA(tag))

  W = Word
  G = Group
  S = Suppress
  O = Optional
  L = Literal
  C = Combine
 
  identifier = W(alphas, alphanums + "_-")
  number = W(nums, nums)
  num_flt = C(number + O(L('.') + number))
  num_sci = C(number + O(L('.') + number) + L('e') + O(L('-') | L('+')) + W(nums))
  gorf = num_sci | num_flt


  wire = G('w' + S('[') + G(number + S(',') + (number | L('f'))) + S(']'))
  gateO = G('g' + S('[') + G(wire + S(',') + number) + S(']'))
  gateI = G('g' + S('[') + G(number+ S(',') + wire) + S(']'))
  thshO = G('th' + S('[') + G(wire + S(',') + number) + S(']'))
  thshI = G('th' + S('[') + G(number + S(',') + wire) + S(']'))

  fluor = G('Fluor' + S('[') + number + S(']'))

  inp = "INPUT" + S(L("(")) + G(number | identifier) + S(")") + S("=") + wire
  out = "OUTPUT" + S("(") + G(number | identifier) + S(")") + S("=") + (fluor | wire)

  inputs = G(S("{") + delimitedList(number, ",") + S("}"))
  outputs= G(S("{") + delimitedList((number | L('f')), ",") + S("}"))
  seesaw = 'seesaw' + S('[') + G(number + S(',') + inputs + S(',') + outputs) + S(']')

  conc = gorf + S(L('*') + L('c'))
  wireconc = 'conc' + S('[') + wire + S(',') + conc + S(']')
  outpconc = 'conc' + S('[') + (gateO | gateI) + S(',') + conc + S(']')
  thshconc = 'conc' + S('[') + (thshO | thshI) + S(',') + conc + S(']')

  # MACROS:
  reporter = 'reporter' + S('[') + G(number + S(',') + number) + S(']')
  inputfanout = 'inputfanout' + S('[') + G(number + S(',') + number + S(',') + inputs) + S(']')
  seesawOR = 'seesawOR' + S('[') + \
      G(number + S(',') + number + S(',') + inputs + S(',') + inputs) + S(']')
  seesawAND = 'seesawAND' + S('[') + \
      G(number + S(',') + number + S(',') + inputs + S(',') + inputs) + S(']')

  macros = reporter | inputfanout | seesawOR | seesawAND
  stmt = G(inp | out | seesaw | wireconc | outpconc | thshconc | macros) + OneOrMore(LineEnd().suppress())

  document = StringStart() + ZeroOrMore(LineEnd().suppress()) + OneOrMore(stmt) + StringEnd()
  document.ignore(pythonStyleComment)

  return document

def parse_seesaw_file(data):
  document = ssw_document_setup()
  return document.parseFile(data).asList()

def parse_seesaw_string(data):
  document = ssw_document_setup()
  return document.parseString(data).asList()
