"""
Augmented Lagrangian and PANOC solvers for nonconvex numerical optimization.
"""
__version__ = "1.0.0a21.dev0"

from .alpaqa import *
from .alpaqa import __c_version__

assert __version__ == __c_version__

from .pyapi import *

# For Sphinx
__all__ = [v for v in dir() if not v.startswith("_") and v != "alpaqa"] + [
    "__c_version__"
]
