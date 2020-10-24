#
# dsdobjects/iupac_utils.py
#   - copy and/or modify together with tests/iupac_utils.py
#
""" Constraint operations using the IUPAC code:
    'A': 'A'
    'C': 'C'
    'G': 'G'
    'T': 'T'
    'R': 'AG'      # purine
    'Y': 'CT'      # pyrimidine
    'S': 'CG'      # strong
    'M': 'AC'
    'W': 'AT'      # weak
    'K': 'GT'
    'V': 'ACG'     # not T
    'H': 'ACT'     # not G
    'D': 'AGT'     # not C
    'B': 'CGT'     # not A
    'N': 'ACGT'
"""

class ConstraintError(Exception):
    pass

wc_complement_dna = {   # ACGT => ACGT
        'A': 'T',       # 1000 => 0001
        'C': 'G',       # 0100 => 0010
        'G': 'C',       # 0010 => 0100
        'T': 'A',       # 0001 => 1000
        'R': 'Y',       # 1010 => 0101
        'Y': 'R',       # 0101 => 1010
        'S': 'B',       # 0110 => 0110
        'M': 'K',       # 1100 => 0011
        'W': 'D',       # 1001 => 1001
        'K': 'N',       # 0011 => 1100
        'V': 'B',       # 1110 => 0111
        'H': 'D',       # 1101 => 1011
        'D': 'N',       # 1011 => 1101
        'B': 'N',       # 0111 => 1110
        'N': 'N'}       # 1111 => 1111

wc_complement_rna = {   # ACGU => ACGU
        'A': 'U',       # 1000 => 0001
        'C': 'G',       # 0100 => 0010
        'G': 'C',       # 0010 => 0100
        'U': 'A',       # 0001 => 1000
        'R': 'Y',       # 1010 => 0101
        'Y': 'R',       # 0101 => 1010
        'S': 'B',       # 0110 => 0110
        'M': 'K',       # 1100 => 0011
        'W': 'D',       # 1001 => 1001
        'K': 'N',       # 0011 => 1100
        'V': 'B',       # 1110 => 0111
        'H': 'D',       # 1101 => 1011
        'D': 'N',       # 1011 => 1101
        'B': 'N',       # 0111 => 1110
        'N': 'N'}       # 1111 => 1111
 
wobble_complement_dna = {# ACGT => ACGT
        'A': 'T',        # 1000 => 0001
        'C': 'G',        # 0100 => 0010
        'G': 'Y',        # 0010 => 0101
        'T': 'R',        # 0001 => 1010
        'R': 'Y',        # 1010 => 0101
        'Y': 'R',        # 0101 => 1010
        'S': 'B',        # 0110 => 0111
        'M': 'K',        # 1100 => 0011
        'W': 'D',        # 1001 => 1011
        'K': 'N',        # 0011 => 1111
        'V': 'B',        # 1110 => 0111
        'H': 'D',        # 1101 => 1011
        'D': 'N',        # 1011 => 1111
        'B': 'N',        # 0111 => 1111
        'N': 'N'}        # 1111 => 1111

wobble_complement_rna = {# ACGU => ACGU
        'A': 'U',       # 1000 => 0001
        'C': 'G',       # 0100 => 0010
        'G': 'Y',       # 0010 => 0101
        'U': 'R',       # 0001 => 1010
        'R': 'Y',       # 1010 => 0101
        'Y': 'R',       # 0101 => 1010
        'S': 'B',       # 0110 => 0111
        'M': 'K',       # 1100 => 0011
        'W': 'D',       # 1001 => 1011
        'K': 'N',       # 0011 => 1111
        'V': 'B',       # 1110 => 0111
        'H': 'D',       # 1101 => 1011
        'D': 'N',       # 1011 => 1111
        'B': 'N',       # 0111 => 1111
        'N': 'N'}       # 1111 => 1111

iupac_bin = {           # ACGT
        'A': 8,         # 1000,
        'C': 4,         # 0100,
        'G': 2,         # 0010,
        'T': 1,         # 0001,
        'U': 1,         # 0001,
        'R': 10,        # 1010,  # purine
        'Y': 5,         # 0101,  # pyrimidine
        'S': 6,         # 0110,
        'M': 12,        # 1100,
        'W': 9,         # 1001,
        'K': 3,         # 0011,
        'V': 14,        # 1110,  # not T/U
        'H': 13,        # 1101,  # not G
        'D': 11,        # 1011,  # not C
        'B': 7,         # 0111,  # not A
        'N': 15}        # 1111,
bin_iupac_dna = ['', 'T', 'G', 'K', 'C', 'Y', 'S', 'B', 'A', 'W', 'R', 'D', 'M', 'H', 'V', 'N']      
bin_iupac_rna = ['', 'U', 'G', 'K', 'C', 'Y', 'S', 'B', 'A', 'W', 'R', 'D', 'M', 'H', 'V', 'N']      


def complement(sequence, material = 'DNA'):
    """ str: complement including wobble base pairs. """
    if material == 'DNA':
        return ''.join([wobble_complement_dna[x] for x in sequence])
    else:
        return ''.join([wobble_complement_rna[x] for x in sequence])

def wc_complement(sequence, material = 'DNA'):
    """ str: Watson-Crick complement. """
    if material == 'DNA':
        return ''.join([wc_complement_dna[x] for x in sequence])
    else:
        return ''.join([wc_complement_rna[x] for x in sequence])

def reverse_complement(sequence, material = 'DNA'):
    """ str: reverse complement including wobble base pairs. """
    if material == 'DNA':
        return ''.join([wobble_complement_dna[x] for x in reversed(sequence)])
    else:
        return ''.join([wobble_complement_rna[x] for x in reversed(sequence)])

def reverse_wc_complement(sequence, material = 'DNA'):
    """ str: reverse Watson-Crick complement. """
    if material == 'DNA':
        return ''.join([wc_complement_dna[x] for x in reversed(sequence)])
    else:
        return ''.join([wc_complement_rna[x] for x in reversed(sequence)])

def add_constraints(seq1, seq2, material = 'DNA'):
    assert len(seq1) == len(seq2)
    if material == 'DNA':
        con = ''.join([bin_iupac_dna[iupac_bin[x] & iupac_bin[y]] for x, y in zip(seq1, seq2)])
    else:
        con = ''.join([bin_iupac_rna[iupac_bin[x] & iupac_bin[y]] for x, y in zip(seq1, seq2)])
    if len(con) < len(seq1):
        raise ConstraintError('Incompatible constraints {seq1} and {seq2}.')

