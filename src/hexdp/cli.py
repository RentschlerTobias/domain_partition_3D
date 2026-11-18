"""Command line interface.

    hexdp --version
    hexdp plugins [KIND]        list registered plugins
    hexdp plugins KIND --describe   load them and print name, doc and parameter schema as JSON

    hexdp run CASE.yaml [--set K=V ...] [--fresh] [--background]   run one case or workflow (resumes)
    hexdp status CASE.yaml                         state of the case's workspace
    hexdp eval CASE.yaml ... [--rerun] [--json F]  run / collect several cases, print a table
    hexdp batch WORKFLOW.yaml --each K=a,b,c [--set K=V] [--workspace D] [--rerun] [--json F]
                                                   run a workflow for several parameter values

    hexdp mcp                                      the MCP server (stdio) for agents
"""

import argparse
import json
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


def _sets(items) -> dict:
    out = {}
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"--set expects NAME=VALUE, got {it!r}")
        k, v = it.split("=", 1)
        out[k] = v
    return out


def _run(args) -> int:
    import json

    from hexdp.pipeline.case import load_case

    overrides = _sets(args.set)
    if args.background:
        from hexdp.pipeline.jobs import start

        print(f"started in background: {start(args.case, overrides)}")
        return 0
    from hexdp.pipeline.runner import run_case

    rep = run_case(load_case(args.case, overrides), fresh=args.fresh)
    keep = {k: rep.get(k) for k in ("case", "parameters", "state", "strategy", "blocks", "reference", "accept",
                                    "accepted", "timings", "error", "time_s") if k in rep}
    print(json.dumps(keep, indent=1, default=str))
    return 0 if rep["state"] == "done" and rep.get("accepted", True) else 1


def _batch(args) -> int:
    import json

    from hexdp.pipeline.evaluate import markdown
    from hexdp.pipeline.workflow import parse_each, run_batch

    runs = parse_each(args.each) if args.each else [{}]
    rows = run_batch(args.workflow, runs, common=_sets(args.set), workspace_root=args.workspace, rerun=args.rerun,
                     status_dir=args.status_dir)
    print(markdown(rows))
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=1, default=str)
    return 0 if all(r.get("state") == "done" and r.get("accept") in ("ok", "-") for r in rows) else 1


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
    r.add_argument("--set", action="append", metavar="NAME=VALUE", help="workflow parameter (repeatable)")
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
    b = sub.add_parser("batch", help="run a workflow for several parameter values")
    b.add_argument("workflow")
    b.add_argument("--each", default=None, metavar="NAME=a,b,c", help="one run per value (or NAME=@file)")
    b.add_argument("--set", action="append", metavar="NAME=VALUE", help="parameter common to all runs")
    b.add_argument("--workspace", default=None, help="put the run workspaces into this directory")
    b.add_argument("--rerun", action="store_true", help="run again even if a finished report exists")
    b.add_argument("--json", default=None, help="write the rows as JSON")
    b.add_argument("--status-dir", default=None, help=argparse.SUPPRESS)   # background validation job
    b.set_defaults(func=_batch)
    xv = sub.add_parser("extract-volume", help="one volume of a multi-volume AlgoHex-format tet mesh (.vtk)")
    xv.add_argument("src")
    xv.add_argument("dst")
    g = xv.add_mutually_exclusive_group(required=True)
    g.add_argument("--volume", type=int, nargs="+", help="volume ids to keep")
    g.add_argument("--role", choices=("outer", "embedded"),
                   help="outer = the volume(s) touching the outer boundary (e.g. box minus body), embedded = the rest")
    xv.set_defaults(func=_extract_volume)
    m = sub.add_parser("mcp", help="run the MCP server (stdio)")
    m.set_defaults(func=lambda a: __import__("hexdp.mcp.server", fromlist=["main"]).main() or 0)
    return p


def _extract_volume(args) -> int:
    from hexdp.io import tetvtk
    m = tetvtk.read(args.src)
    vols = args.volume or tetvtk.volume_roles(m)[args.role]
    if not vols:
        print(f"no {args.role} volume (volumes: {m.volumes()})", file=sys.stderr)
        return 1
    sub = tetvtk.extract_volume(m, vols)
    tetvtk.write(sub, args.dst)
    print(json.dumps({"volumes": vols, "tets": len(sub.tets), "boundary_triangles": len(sub.tris),
                      "face_ids": sorted(set(sub.tri_ids.tolist())), "dst": str(args.dst)}))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
