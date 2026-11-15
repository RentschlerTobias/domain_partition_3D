"""Command line interface.

    hexdp --version
    hexdp plugins [KIND]        list registered plugins
    hexdp plugins KIND --describe   load them and print name, doc and parameter schema as JSON

    hexdp run CASE.yaml [--fresh] [--background]   run one case (resumes its workspace)
    hexdp status CASE.yaml                         state of the case's workspace
    hexdp eval CASE.yaml ... [--rerun] [--json F]  run / collect several cases, print a table

    hexdp mcp                                      the MCP server (stdio) for agents
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


def _run(args) -> int:
    import json

    from hexdp.pipeline.case import load_case

    if args.background:
        from hexdp.pipeline.jobs import start

        print(f"started in background: {start(args.case)}")
        return 0
    from hexdp.pipeline.runner import run_case

    rep = run_case(load_case(args.case), fresh=args.fresh)
    keep = {k: rep.get(k) for k in ("case", "state", "strategy", "blocks", "reference", "timings", "error", "time_s")}
    print(json.dumps(keep, indent=1, default=str))
    return 0 if rep["state"] == "done" else 1


def _status(args) -> int:
    import json

    from hexdp.pipeline.case import load_case
    from hexdp.pipeline.runner import status

    case = load_case(args.case)
    print(json.dumps({"workspace": str(case.workspace), **status(case.workspace)}, indent=1, default=str))
    return 0


def _eval(args) -> int:
    import json

    from hexdp.pipeline.evaluate import evaluate, markdown

    rows = evaluate(args.cases, rerun=args.rerun)
    print(markdown(rows))
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=1, default=str)
    return 0 if all(r.get("state") == "done" and not str(r.get("reference", "")).startswith("FAIL") for r in rows) else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="hexdp", description="Hexahedral block structures and domain partition.")
    p.add_argument("--version", action="version", version=f"hexdp {__version__}")
    sub = p.add_subparsers(dest="command")
    pl = sub.add_parser("plugins", help="list registered plugins")
    pl.add_argument("kind", nargs="?", choices=registry.KINDS)
    pl.add_argument("--describe", action="store_true", help="load and describe (JSON)")
    pl.set_defaults(func=_plugins)
    r = sub.add_parser("run", help="run one case file")
    r.add_argument("case")
    r.add_argument("--fresh", action="store_true", help="refuse a non-empty workspace")
    r.add_argument("--background", action="store_true", help="start detached, return at once")
    r.set_defaults(func=_run)
    s_ = sub.add_parser("status", help="state of a case's workspace")
    s_.add_argument("case")
    s_.set_defaults(func=_status)
    e = sub.add_parser("eval", help="run or collect several cases and summarise")
    e.add_argument("cases", nargs="+")
    e.add_argument("--rerun", action="store_true", help="run again even if a finished report exists")
    e.add_argument("--json", default=None, help="write the rows as JSON")
    e.set_defaults(func=_eval)
    m = sub.add_parser("mcp", help="run the MCP server (stdio)")
    m.set_defaults(func=lambda a: __import__("hexdp.mcp.server", fromlist=["main"]).main() or 0)
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
