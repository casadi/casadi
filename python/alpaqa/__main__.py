import argparse
from .cache import get_alpaqa_cache_dir, interactive_clean
from . import __version__


def main():
    # alpaqa
    parser = argparse.ArgumentParser(description="alpaqa cache operations",
                                     allow_abbrev=False)
    parser.add_argument('--version',
                        '-v',
                        action='store_true',
                        help='print the version number and exit')
    subparsers = parser.add_subparsers(dest="command")
    # alpaqa cache
    cacheparser = subparsers.add_parser("cache")
    cachesubparsers = cacheparser.add_subparsers(dest="cachecommand")
    # alpaqa cache clean
    cleanparser = cachesubparsers.add_parser(
        "clean", description="clean up all cached files")
    cleanparser.add_argument(
        '--force',
        '-f',
        action='store_true',
        help='do not ask for confirmation before deleting')
    cleanparser.add_argument(
        '--dry',
        '-n',
        action='store_true',
        help='perform a dry run without actually deleting the cache')
    # alpaqa cache print
    cachesubparsers.add_parser("print",
                               description="print the cache path and exit")

    args = parser.parse_args()

    if args.version:
        print(__version__)
    elif args.command == 'cache':
        if args.cachecommand == 'clean':
            interactive_clean(args)
        elif args.cachecommand == 'print':
            print(get_alpaqa_cache_dir())
            exit(0)


if __name__ == '__main__':
    main()
