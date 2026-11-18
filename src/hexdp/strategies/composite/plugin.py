"""Strategy plugin ``composite``: a block structure assembled from building blocks (plugin kind ``blocks``) that a
workflow lists as ``steps``:

    strategy:
      - name: composite
        params:
          steps:
            - {id: rotor, use: dtoo_rotor, params: {fillet: true}}          # inputs: the case geometry
            - {id: hub_surface, use: band_layout, from: [rotor], params: {h: 0.04}}
            - {id: hub, use: revolve_fill, from: [hub_surface], params: {r_bore: 0.15}}
            - {id: blades, use: section_sweep, from: [hub_surface, rotor], params: {h: 0.04}}

Every step gets the results of the steps named in ``from`` (and ``geometry``, the case geometry) and returns a
result; the results that are parts (hexdp.blocks.Part) are composed into one conforming block structure (coincident
fine nodes merged), boundary faces named by their parts.
"""

from __future__ import annotations

import json
import time

from hexdp import registry
from hexdp.blocks import Part, compose
from hexdp.contracts import ParamSpec, Score, resolve_params


class CompositeStrategy:
    """Block structure from building blocks: input, surface layouts, fills and sweeps per sub-domain, composed at
    conforming interfaces. The steps come from the workflow (params.steps)."""

    name = "composite"
    params = {"steps": ParamSpec([], "list of {id, use, from, params}: building blocks in order")}

    def applicable(self, geo, features):
        return Score(0.3, ["composes building blocks as the workflow lists them (needs params.steps)"])

    def run(self, geo, features, params, ctx):
        steps = params["steps"]
        if not steps:
            raise ValueError("composite: no steps given (params.steps)")
        results, parts, log = {"geometry": geo}, [], []
        for st in steps:
            sid, use = st["id"], st["use"]
            plugin = registry.load("blocks", use)
            src = st.get("from", [])
            src = [src] if isinstance(src, str) else list(src)
            missing = [s for s in src if s not in results]
            if missing:
                raise ValueError(f"step {sid}: unknown input(s) {missing}; known: {sorted(results)}")
            inputs = {"geometry": geo, **{s: results[s] for s in src}}
            t = time.time()
            out = plugin.run(inputs, resolve_params(plugin.params, st.get("params")), ctx)
            results[sid] = out
            if isinstance(out, Part):
                parts.append(out)
            log.append({"id": sid, "use": use, "time_s": round(time.time() - t, 2),
                        **({"part": out.info} if isinstance(out, Part) else {})})
            ctx.decide("composite", f"{sid} ({use})", [json.dumps(log[-1], default=float)])
        if not parts:
            raise ValueError("composite: no step produced a part")
        bs, info = compose(parts)
        (ctx.workspace / "composite.json").write_text(json.dumps({"steps": log, **info}, indent=1, default=float))
        return bs
