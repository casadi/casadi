import os
import typing


def _is_truthy(s: typing.Optional[str]):
    if s is None:
        return False
    return not s.lower() in ("", "false", "no", "off", "0")


if not typing.TYPE_CHECKING and _is_truthy(os.getenv("ALPAQA_PYTHON_DEBUG")):
    from . import _alpaqa_d
    from ._alpaqa_d import *
    from ._alpaqa_d.float64 import *
    from ._alpaqa_d import __version__ as __c_version__
else:
    from . import _alpaqa
    from ._alpaqa import *
    from ._alpaqa.float64 import *
    from ._alpaqa import __version__ as __c_version__

del _is_truthy, typing, os
