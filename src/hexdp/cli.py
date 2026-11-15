"""Command line interface.

    hexdp --version
    hexdp plugins [KIND]        list registered plugins
    hexdp plugins KIND --describe   load them and print name, doc and parameter schema as JSON

The commands ``run`` (one case YAML), ``eval`` (a whole dataset) and ``mcp``
(the MCP server) are added in steps 1.5 and 1.6 of the plan.
"""

import argparse
import sys

from hexdp import __version__, registry


def _plugins(args: argparse.Namespace) -> int:
    if args.describe:
        import json

        from hexdp.contracts import describe

        kinds = [args.kind] if args.kind else list(registry.KINDS)
        out = {k: [describe(registry.load(k, n)) for n in registry.names(k)] for k in kinds}
        print(json.dumps(out, indent=1))
        return 0
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
    pl.add_argument("--describe", action="store_true", help="load and describe (JSON)")
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
