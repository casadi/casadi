import os
from os.path import join, expanduser
import tempfile
import shutil


def get_cache_dir():
    """
    Get the temporary directory used for caching compiled problems.

    If the environment variable ``ALPAQA_CACHE_DIR`` exists, its value is
    returned.
    Otherwise, if ``~/.cache`` exists, then it is returned.
    Otherwise, ``tempfile.gettempdir()`` is returned.
    """
    cache_dir = os.getenv("ALPAQA_CACHE_DIR")
    if cache_dir is None:
        home_cache_dir = expanduser("~/.cache")
        if os.path.isdir(home_cache_dir):
            cache_dir = home_cache_dir
        else:
            cache_dir = tempfile.gettempdir()
    return cache_dir


def get_alpaqa_cache_dir():
    return join(get_cache_dir(), 'alpaqa', 'cache')


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
            confirmed = input(msg).lower() == 'y'
        if not confirmed:
            print("Aborted.")
            exit(1)
        print("Removing", alpaqa_cache_dir)
        clean(alpaqa_cache_dir)
        exit(0)
