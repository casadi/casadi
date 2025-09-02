"""
Augmented Lagrangian and PANOC solvers for nonconvex numerical optimization.
"""

__version__ = "1.1.0a1"

import contextlib
from .alpaqa import *
from .alpaqa import __c_version__

assert __version__ == __c_version__

with contextlib.suppress(ModuleNotFoundError):  # Don't fail if CasADi is unavailable
    from .pyapi import *

# For Sphinx
__all__ = [v for v in dir() if not v.startswith("_") and v != "alpaqa"] + [
    "__c_version__"
]
