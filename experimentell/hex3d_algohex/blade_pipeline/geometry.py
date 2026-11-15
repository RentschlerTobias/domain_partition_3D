"""Moved to `hexdp.strategies.blade_ring.geometry` (src/hexdp/strategies/blade_ring/geometry.py).

This shim keeps the old imports (`import geometry`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/geometry.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.blade_ring.geometry", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.blade_ring.geometry as _moved

    sys.modules[__name__] = _moved
