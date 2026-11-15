"""The stage chain: strategy -> ops -> discretiser -> filler -> exporters, with decisions logged.

This is the minimal composition the contracts promise; the runner with case files, jobs and
timeouts (step 1.5) builds on it. Strategy choice: the given name, or with ``"auto"`` the applicable
strategy with the highest score.
"""

from __future__ import annotations

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


def run_chain(geo: Geometry, ctx: Context, *, detectors=(), strategies, ops=(), discretiser, filler,
              exporters=()) -> ChainResult:
    feats = detect(geo, list(detectors), ctx)
    st, scores = choose_strategy(geo, feats, list(strategies), ctx)
    bs = st.plugin.run(geo, feats, resolve_params(st.plugin.params, st.params), ctx)
    for op in ops:
        sc = op.plugin.applicable(bs, geo)
        if not sc.ok:
            ctx.decide("op", f"skip {op.plugin.name}", sc.reasons)
            continue
        bs = op.plugin.apply(bs, geo, resolve_params(op.plugin.params, op.params), ctx)
        ctx.decide("op", op.plugin.name, sc.reasons + [f"{bs.n_blocks} blocks after"])
    bs = discretiser.plugin.apply(bs, geo, resolve_params(discretiser.plugin.params, discretiser.params), ctx)
    bs = filler.plugin.fill(bs, geo, resolve_params(filler.plugin.params, filler.params), ctx)
    for ex in exporters:
        written = ex.plugin.write(bs, ctx.path(f"{ex.plugin.name}"), resolve_params(ex.plugin.params, ex.params), ctx)
        ctx.decide("export", ex.plugin.name, [str(p) for p in written])
    ctx.save_decisions()
    return ChainResult(bs, feats, scores, st.plugin.name)
