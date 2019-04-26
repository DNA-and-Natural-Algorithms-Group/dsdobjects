#
# dsdobjects.parser.pil_extended_format
#   - copy and/or modify together with tests/test_pil_parser.py
#
# Written by Stefan Badelt (badelt@caltech.edu)
#
# Distributed under the MIT License, use at your own risk.
#

from pyparsing import (Word, Literal, Group, Suppress, Optional, ZeroOrMore,
        Combine, White, OneOrMore, alphas, alphanums, nums, delimitedList,
        StringStart, StringEnd, Forward, LineEnd, pythonStyleComment,
        ParseElementEnhance)

class PilFormatError(Exception):
    pass

def pil_extended_setup():
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
    C = Combine
    L = Literal

    identifier = W(alphanums + "_-")
    number = W(nums, nums)

    num_flt = C(number + O(L('.') + number))
    num_sci = C(number + O(L('.') + number) + L('e') + O(L('-') | L('+')) + W(nums))
    gorf = num_sci | num_flt
    ginf = gorf | L('inf')

    domain = C(identifier + O('*'))
    constraint = W(alphas)
    assign = L("=") | L(":")
    dotbracket = W("(.)+ ")
    dlength = number | L('short') | L('long')

    dl_domain = G(T(S("length") + domain + S(assign) + dlength + OneOrMore(LineEnd().suppress()), 'dl-domain')) \
              | G(T(S("domain") + domain + S(assign) + dlength + OneOrMore(LineEnd().suppress()), 'dl-domain')) \
              | G(T(S("sequence") + domain + S(assign) + dlength + OneOrMore(LineEnd().suppress()), 'dl-domain'))

    sl_domain = G(T(S("sequence") + domain + S(assign) + constraint + O(S(assign) + number) + OneOrMore(LineEnd().suppress()), 'sl-domain'))

    # strand and sup-sequence are the same thing ...
    comp_domain = G(T(S("sup-sequence") + identifier + S(assign) \
            + G(OneOrMore(domain)) + O(S(assign) + number) \
            + OneOrMore(LineEnd().suppress()), 'composite-domain'))
    strand = G(T(S("strand") + identifier + S(assign) \
            + G(OneOrMore(domain)) + O(S(assign) + number) \
            + OneOrMore(LineEnd().suppress()), 'composite-domain'))

    strandcomplex = G(T(S("complex") + identifier + S(assign) + O(LineEnd().suppress()) \
                    + G(OneOrMore(domain)) + O(LineEnd().suppress()) \
                    + dotbracket + OneOrMore(LineEnd().suppress()), 'strand-complex')) \
                  | G(T(S("structure") + identifier + S(assign) \
                    + G(OneOrMore(domain | S('+'))) + S(assign) \
                    + dotbracket + OneOrMore(LineEnd().suppress()), 'strand-complex'))

    species = delimitedList(identifier, '+')
    cunit = L('M') | L('mM') | L('uM') | L('nM') | L('pM') 
    tunit = L('s') | L('m') | L('h')
    runit = C(ZeroOrMore('/' + cunit) + L('/') + tunit)
    infobox = S('[') + G(O(identifier + S(assign))) + G(gorf + O(S(L('+/-')) + ginf)) + G(runit) + S(']')

    reaction = G(T(S("kinetic") + G(O(infobox)) + G(species) + S('->') + G(species) + OneOrMore(LineEnd().suppress()), 'reaction')) \
             | G(T(S("reaction") + G(O(infobox)) + G(species) + S('->') + G(species) + OneOrMore(LineEnd().suppress()), 'reaction'))

    restingset = G(T(S("state") + identifier + S("=") + S('[') + G(delimitedList(identifier)) + S(']') + OneOrMore(LineEnd().suppress()), 'resting-macrostate')) \
               | G(T(S("macrostate") + identifier + S("=") + S('[') + G(delimitedList(identifier)) + S(']') + OneOrMore(LineEnd().suppress()), 'resting-macrostate'))

    # kernel notation
    sense = Combine(identifier + O(L("^")) + O(L("*")))

    pattern = Forward()
    innerloop = pattern | S(White())
    loop = (Combine(sense + S("(")) + G(O(innerloop)) + S(")"))
    pattern << OneOrMore(loop | L("+") | sense)

    conc = G( S('@') + (L('initial')  | L('i')) + gorf + cunit) \
         | G( S('@') + (L('constant') | L('c')) + gorf + cunit)

    cplx = G(T(identifier + S("=") + OneOrMore(G(pattern)) + O(conc) +
        OneOrMore(LineEnd().suppress()), 'kernel-complex'))

    stmt = sl_domain | dl_domain | comp_domain | strand | strandcomplex | reaction | cplx | restingset

    document = StringStart() + ZeroOrMore(LineEnd().suppress()) + \
        OneOrMore(stmt) + StringEnd()
    document.ignore(pythonStyleComment)

    return document

def parse_pil_file(data):
    document = pil_extended_setup()
    return document.parseFile(data).asList()

def parse_pil_string(data):
    document = pil_extended_setup()
    return document.parseString(data).asList()

