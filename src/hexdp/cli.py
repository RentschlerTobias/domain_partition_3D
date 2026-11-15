"""Command line interface.

    hexdp --version
    hexdp plugins [KIND]        list registered plugins

The commands ``run`` (one case YAML), ``eval`` (a whole dataset) and ``mcp``
(the MCP server) are added in steps 1.5 and 1.6 of the plan.
"""

import argparse
import sys

from hexdp import __version__, registry


def _plugins(args: argparse.Namespace) -> int:
    kinds = [args.kind] if args.kind else list(registry.KINDS)
    for kind in kinds:
        found = registry.entries(kind)
        print(f"{kind}:" + ("" if found else " (none)"))
        for e in found:
            origin = f"  [{e.distribution}]" if e.distribution else ""
            print(f"  {e.name:20s} {e.target}{origin}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="hexdp", description="Hexahedral block structures and domain partition.")
    p.add_argument("--version", action="version", version=f"hexdp {__version__}")
    sub = p.add_subparsers(dest="command")
    pl = sub.add_parser("plugins", help="list registered plugins")
    pl.add_argument("kind", nargs="?", choices=registry.KINDS)
    pl.set_defaults(func=_plugins)
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
