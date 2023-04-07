"""
Augmented Lagrangian and PANOC solvers for nonconvex numerical optimization.
"""
__version__ = "1.0.0a7"

from .alpaqa import *
from .alpaqa import __c_version__

assert __version__ == __c_version__

from .pyapi import *
