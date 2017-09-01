#
#
# Copyright (c) 2017 Caltech. All rights reserved.
# Written by Stefan Badelt (badelt@caltech.edu).
#
#
# Parser module for pil files in kernel notation.
#

from pyparsing import (Word, Literal, Group, Suppress, Optional, ZeroOrMore,
        Combine, White, OneOrMore, alphas, alphanums, nums, StringStart,
        StringEnd, Forward, LineEnd, pythonStyleComment, ParseElementEnhance)


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
    gorf = C(W(nums) + O((L('.') + W(nums)) | (L('e') + O('-') + W(nums)))) # {:g} {:f}
    domain = G(T(S("length") + identifier + S("=") + number +
                 OneOrMore(LineEnd().suppress()), 'domain'))

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

    stmt = domain | cplx

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

