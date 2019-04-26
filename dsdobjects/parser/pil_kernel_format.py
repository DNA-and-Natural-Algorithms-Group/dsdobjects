#
# dsdobjects.parser.pil_kernel_format
#   - copy and/or modify together with tests/test_kernel_parser.py
#
# Written by Stefan Badelt (badelt@caltech.edu)
#
# Distributed under the MIT License, use at your own risk.
#

from pyparsing import (Word, Literal, Group, Suppress, Optional, ZeroOrMore,
        Combine, White, OneOrMore, alphas, alphanums, nums, delimitedList,
        StringStart, StringEnd, Forward, LineEnd, pythonStyleComment,
        ParseElementEnhance)


def pil_kernel_setup():
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

    identifier = W(alphas, alphanums + "_-") # forbid names starting with digits
    number = W(nums, nums)

    num_flt = C(number + O(L('.') + number))
    num_sci = C(number + O(L('.') + number) + L('e') + O(L('-') | L('+')) + W(nums))
    gorf = num_sci | num_flt

    dlength = number | L('short') | L('long')

    domain = G(T(S("length") + identifier + S("=") + dlength +
                 OneOrMore(LineEnd().suppress()), 'dl-domain'))

    # NOTE: exchange the comment for asense if you want to allow input in form
    # of "x( ... y)", but also double-check if that really works...
    sense = Combine(identifier + O(L("^")) + O(L("*")))
    sbreak = L("+")

    pattern = Forward()
    innerloop = S(White()) + pattern + S(White()) | G(S(White()))
    loop = (Combine(sense + S("(")) + innerloop + S(")"))
    pattern << G(OneOrMore(loop | sbreak | sense))

    unit = L('M') | L('mM') | L('uM') | L('nM') | L('pM')
    conc = G( S('@') + L('initial') + gorf + unit) \
         | G( S('@') + L('constant') + gorf + unit)

    cplx = G(T(identifier + S("=") + OneOrMore(pattern) +
               O(conc) + OneOrMore(LineEnd().suppress()), 'complex'))

    restingset = G(T(S("state") + identifier + S("=") + S('[') + G(delimitedList(identifier)) + S(']') + OneOrMore(LineEnd().suppress()), 'resting-macrostate')) \
               | G(T(S("macrostate") + identifier + S("=") + S('[') + G(delimitedList(identifier)) + S(']') + OneOrMore(LineEnd().suppress()), 'resting-macrostate'))

    species = delimitedList(identifier, '+')
    units = W("/M/s")
    infobox = S('[') + G(O(identifier + S('='))) + G(gorf) + G(units) + S(']')
    reaction = G(T(S("reaction") + G(O(infobox)) + G(species) + S('->') + G(species) + OneOrMore(LineEnd().suppress()), 'reaction'))

    stmt = domain | cplx | reaction | restingset

    document = StringStart() + ZeroOrMore(LineEnd().suppress()) + \
        OneOrMore(stmt) + StringEnd()
    document.ignore(pythonStyleComment)

    return document


def parse_kernel_file(data):
    document = pil_kernel_setup()
    return document.parseFile(data).asList()


def parse_kernel_string(data):
    document = pil_kernel_setup()
    return document.parseString(data).asList()

