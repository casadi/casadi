import os
from os.path import join, expanduser
import tempfile
import shutil
import sys
import typing


def _is_truthy(s: typing.Optional[str]):
    if s is None:
        return False
    return not s.lower() in ("", "false", "no", "off", "0")


def get_cache_dir():
    """
    Get the temporary directory used for caching compiled problems.

    If the environment variable ``ALPAQA_CACHE_DIR`` exists, its value is
    returned.
    Otherwise, if a virtual environment is active and the environment variable
    ``ALPAQA_GLOBAL_CACHE`` is unset (or set to a false value),
    ``{sys.prefix}/cache`` is returned.
    Otherwise, if ``~/.cache`` exists, then it is returned.
    Otherwise, the result of ``tempfile.gettempdir()`` is returned.
    """
    # 1. Custom directory
    cache_dir = os.getenv("ALPAQA_CACHE_DIR")
    if cache_dir is not None:
        return cache_dir
    # 2. Virtual environment cache directory
    in_venv = sys.prefix != sys.base_prefix
    if in_venv and not _is_truthy(os.getenv("ALPAQA_GLOBAL_CACHE")):
        return join(sys.prefix, "cache")
    # 3. Home cache directory
    home_cache_dir = expanduser("~/.cache")
    if os.path.isdir(home_cache_dir):
        return home_cache_dir
    # 4. System temporary directory
    return tempfile.gettempdir()


def get_alpaqa_cache_dir():
    return join(get_cache_dir(), "alpaqa", "cache")


def clean(alpaqa_cache_dir=None):
    if alpaqa_cache_dir is None:
        alpaqa_cache_dir = get_alpaqa_cache_dir()
    shutil.rmtree(alpaqa_cache_dir, ignore_errors=True)


def interactive_clean(args):
    alpaqa_cache_dir = get_alpaqa_cache_dir()

    if args.dry:
        print(f"Would remove: {alpaqa_cache_dir}")
        exit(0)
    else:
        confirmed = args.force
        msg = f"Are you sure you want to remove {alpaqa_cache_dir} [y/N]? "
        if not confirmed:
            confirmed = input(msg).lower() == "y"
        if not confirmed:
            print("Aborted.")
            exit(1)
        print("Removing", alpaqa_cache_dir)
        clean(alpaqa_cache_dir)
        exit(0)
