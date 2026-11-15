"""Moved to `hexdp.strategies.blade_ring.snap` (src/hexdp/strategies/blade_ring/snap.py).

This shim keeps the old imports (`import snap`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/snap.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.blade_ring.snap", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.blade_ring.snap as _moved

    sys.modules[__name__] = _moved
