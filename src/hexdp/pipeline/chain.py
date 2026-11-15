"""The stage chain: strategy -> ops -> discretiser -> filler -> exporters, with decisions logged.

This is the minimal composition the contracts promise; the runner with case files, jobs and
timeouts (step 1.5) builds on it. Strategy choice: the given name, or with ``"auto"`` the applicable
strategy with the highest score.
"""

from __future__ import annotations

import time

from dataclasses import dataclass, field

from hexdp.contracts import Context, Feature, Score, resolve_params
from hexdp.core.blockstructure import BlockStructure
from hexdp.core.geometry import Geometry


@dataclass
class Stage:
    plugin: object
    params: dict = field(default_factory=dict)


@dataclass
class ChainResult:
    blocks: BlockStructure
    features: list[Feature]
    scores: dict[str, Score]
    strategy: str
    timings: dict[str, float] = field(default_factory=dict)
    outputs: dict[str, list[str]] = field(default_factory=dict)


def detect(geo: Geometry, detectors: list[Stage], ctx: Context) -> list[Feature]:
    feats: list[Feature] = []
    for st in detectors:
        found = st.plugin.detect(geo, resolve_params(st.plugin.params, st.params), ctx)
        ctx.decide("detect", st.plugin.name, [f"{f.kind} ({f.confidence:.2f}): {'; '.join(f.reasons)}" for f in found]
                   or ["nothing found"])
        feats += found
    return feats


def choose_strategy(geo: Geometry, feats: list[Feature], strategies: list[Stage], ctx: Context):
    scores = {st.plugin.name: st.plugin.applicable(geo, feats) for st in strategies}
    ok = [st for st in strategies if scores[st.plugin.name].ok]
    if not ok:
        why = {n: s.reasons + [f"needs {m}" for m in s.missing] for n, s in scores.items()}
        raise RuntimeError(f"no applicable strategy: {why}")
    best = max(ok, key=lambda st: scores[st.plugin.name].value)
    ctx.decide("strategy", best.plugin.name, scores[best.plugin.name].reasons,
               {n: s.value for n, s in scores.items()})
    return best, scores


def run_chain(geo: Geometry, ctx: Context, *, detectors=(), strategies, ops=(), discretiser=None, filler=None,
              exporters=()) -> ChainResult:
    """Run the stages. Discretiser and filler are optional: some strategies deliver a filled structure."""
    timings: dict[str, float] = {}
    outputs: dict[str, list[str]] = {}

    def timed(name, fn, *a):
        t = time.time()
        try:
            return fn(*a)
        finally:
            timings[name] = round(time.time() - t, 3)

    feats = timed("detect", detect, geo, list(detectors), ctx)
    st, scores = choose_strategy(geo, feats, list(strategies), ctx)
    bs = timed(f"strategy:{st.plugin.name}", st.plugin.run, geo, feats, resolve_params(st.plugin.params, st.params), ctx)
    for op in ops:
        sc = op.plugin.applicable(bs, geo)
        if not sc.ok:
            ctx.decide("op", f"skip {op.plugin.name}", sc.reasons + [f"needs {m}" for m in sc.missing])
            continue
        n0 = bs.n_blocks
        bs = timed(f"op:{op.plugin.name}", op.plugin.apply, bs, geo, resolve_params(op.plugin.params, op.params), ctx)
        ctx.decide("op", op.plugin.name, sc.reasons + [f"{n0} -> {bs.n_blocks} blocks"])
    if discretiser is not None:
        bs = timed(f"discretise:{discretiser.plugin.name}", discretiser.plugin.apply, bs, geo,
                   resolve_params(discretiser.plugin.params, discretiser.params), ctx)
    if filler is not None:
        bs = timed(f"fill:{filler.plugin.name}", filler.plugin.fill, bs, geo,
                   resolve_params(filler.plugin.params, filler.params), ctx)
    for ex in exporters:
        written = timed(f"export:{ex.plugin.name}", ex.plugin.write, bs, ctx.path("out", ex.plugin.name),
                        resolve_params(ex.plugin.params, ex.params), ctx)
        outputs[ex.plugin.name] = [str(p) for p in written]
        ctx.decide("export", ex.plugin.name, outputs[ex.plugin.name])
    ctx.save_decisions()
    return ChainResult(bs, feats, scores, st.plugin.name, timings, outputs)
