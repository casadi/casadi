"""
Augmented Lagrangian and PANOC solvers for nonconvex numerical optimization.
"""
__version__ = '1.0.0a0'

import os
import typing

def _is_truthy(s: typing.Optional[str]):
    if s is None: return False
    return not s.lower() in ('', 'false', 'no', 'off', '0')

if not typing.TYPE_CHECKING and _is_truthy(os.getenv('ALPAQA_PYTHON_DEBUG')):
    from alpaqa._alpaqa_d import *
    from alpaqa._alpaqa_d.float64 import *
    from alpaqa._alpaqa_d import __version__ as __c_version__
else:
    from alpaqa._alpaqa import *
    from alpaqa._alpaqa.float64 import *
    from alpaqa._alpaqa import __version__ as __c_version__
assert __version__ == __c_version__

del _is_truthy, typing, os