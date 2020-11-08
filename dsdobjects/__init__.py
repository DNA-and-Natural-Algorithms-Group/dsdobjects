#
# dsdobjects
#
#  Contributions:
#  This package contains adapted code from various related Python packages
#  coded in the "DNA and Natural Algorithms Group", Caltech:
#   * "DNAObjecs" coded by Joseph Berleant and Joseph Schaeffer 
#   * "Peppercornenumerator" coded by Kathrik Sarma, Casey Grun and Erik Winfree
#   * "Nuskell" coded by Seung Woo Shin
#
__version__='0.8'

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

from .dsdparser import *
from .singleton import (SingletonError, 
                        clear_singletons,
                        show_singletons)
from .base_classes import (ObjectInitError, 
                           DomainS, 
                           ComplexS,
                           MacrostateS,
                           ReactionS) 
from .iupac_utils import ConstraintError
from .complex_utils import SecondaryStructureError
from .objectio import (read_pil, read_pil_line)

# Deprecated since v0.8, 
from .core.deprecated import clear_memory, DSDObjectsError, DSDDuplicationError
