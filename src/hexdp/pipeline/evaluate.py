"""Evaluate a set of cases: run each (or reuse a finished report) and summarise the results."""

from __future__ import annotations

import json
from pathlib import Path

from hexdp.pipeline.case import CaseError, load_case
from hexdp.pipeline.runner import run_case

COLUMNS = ("case", "state", "strategy", "blocks", "min_sj", "issues", "time_s", "reference")


def row(report: dict) -> dict:
    b = report.get("blocks") or {}
    ref = report.get("reference") or {}
    return {
        "case": report.get("case"),
        "state": report.get("state"),
        "strategy": report.get("strategy"),
        "blocks": b.get("blocks"),
        "min_sj": None if b.get("min_sj_fine") is None else round(b["min_sj_fine"], 4),
        "issues": ",".join(f"{k}:{v}" for k, v in (b.get("issues") or {}).items()) or "-",
        "time_s": report.get("time_s"),
        "reference": "-" if not ref else ("ok" if all(v["ok"] for v in ref.values()) else
                                          "FAIL " + ",".join(k for k, v in ref.items() if not v["ok"])),
        "error": report.get("error"),
    }


def evaluate(case_files, rerun: bool = False) -> list[dict]:
    rows = []
    for f in case_files:
        try:
            case = load_case(f)
        except CaseError as exc:
            rows.append({"case": Path(f).stem, "state": "invalid", "error": str(exc)})
            continue
        rep_file = case.workspace / "report.json"
        if not rerun and rep_file.exists() and json.loads(rep_file.read_text()).get("state") == "done":
            report = json.loads(rep_file.read_text())
        else:
            report = run_case(case)
        rows.append(row(report))
    return rows


def markdown(rows) -> str:
    out = ["| " + " | ".join(COLUMNS) + " |", "|" + "---|" * len(COLUMNS)]
    for r in rows:
        out.append("| " + " | ".join("" if r.get(c) is None else str(r.get(c)) for c in COLUMNS) + " |")
    errs = [f"- {r['case']}: {r['error']}" for r in rows if r.get("error")]
    return "\n".join(out + ([""] + errs if errs else []))
