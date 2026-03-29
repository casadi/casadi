#!/usr/bin/env python3
"""Generate a broad stub surface and a pyright access check from a live module.

The generated output is intentionally wide: it covers every exported module name
plus the Python protocol methods that make CasADi's operator-heavy classes usable.
"""

from __future__ import annotations

import argparse
import difflib
import importlib
from pathlib import Path
from typing import Any


# These dunder methods are part of the Python-facing API for CasADi objects even
# though they don't look like normal "public" names.
SPECIAL_METHODS = frozenset(
    {
        "__abs__",
        "__add__",
        "__array__",
        "__bool__",
        "__call__",
        "__contains__",
        "__eq__",
        "__float__",
        "__ge__",
        "__getitem__",
        "__gt__",
        "__index__",
        "__int__",
        "__iter__",
        "__le__",
        "__len__",
        "__lt__",
        "__matmul__",
        "__mul__",
        "__ne__",
        "__neg__",
        "__next__",
        "__pos__",
        "__pow__",
        "__radd__",
        "__rmatmul__",
        "__rmul__",
        "__rpow__",
        "__rsub__",
        "__rtruediv__",
        "__setitem__",
        "__sub__",
        "__truediv__",
    }
)


def exported_names(obj: Any) -> list[str]:
    return sorted(
        name
        for name in dir(obj)
        if not name.startswith("_") or name in SPECIAL_METHODS
    )


def render_class(name: str, cls: type[Any]) -> list[str]:
    lines = [f"class {name}:"]
    members = exported_names(cls)
    if not members:
        lines.append("    ...")
        return lines

    for member in members:
        try:
            attr = getattr(cls, member)
        except Exception:
            attr = None

        if callable(attr):
            lines.append(f"    def {member}(self, *args: Any, **kwargs: Any) -> Any: ...")
        else:
            lines.append(f"    {member}: Any")

    return lines


def render_stub(module_name: str) -> str:
    module = importlib.import_module(module_name)
    lines = [
        "from __future__ import annotations",
        "",
        "from typing import Any",
        "",
    ]

    for name in exported_names(module):
        obj = getattr(module, name)
        if isinstance(obj, type):
            lines.extend(render_class(name, obj))
        elif callable(obj):
            lines.append(f"def {name}(*args: Any, **kwargs: Any) -> Any: ...")
        else:
            lines.append(f"{name}: Any")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def render_access_check(module_name: str) -> str:
    module = importlib.import_module(module_name)
    lines = [
        "from __future__ import annotations",
        "",
        "from typing import Any",
        "",
        f"import {module_name} as ca",
        "",
        "def touch(value: Any) -> None: ...",
        "",
    ]

    for name in exported_names(module):
        lines.append(f"touch(ca.{name})")
        obj = getattr(module, name)
        if isinstance(obj, type):
            for member in exported_names(obj):
                lines.append(f"touch(ca.{name}.{member})")

    return "\n".join(lines).rstrip() + "\n"


def update_file(content: str, path: Path, check: bool) -> int:
    if check:
        current = path.read_text() if path.exists() else ""
        if current == content:
            return 0

        diff = difflib.unified_diff(
            current.splitlines(),
            content.splitlines(),
            fromfile=str(path),
            tofile=f"{path} (generated)",
            lineterm="",
        )
        print("\n".join(diff))
        return 1

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["stub", "access-check"])
    parser.add_argument("--module", default="casadi")
    parser.add_argument("--output", required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    if args.mode == "stub":
        content = render_stub(args.module)
    else:
        content = render_access_check(args.module)

    return update_file(content, Path(args.output), args.check)


if __name__ == "__main__":
    raise SystemExit(main())
