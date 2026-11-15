"""Moved to `hexdp.strategies.blade_ring.extract_blade` (src/hexdp/strategies/blade_ring/extract_blade.py).

This shim keeps the old imports (`import extract_blade`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/extract_blade.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.blade_ring.extract_blade", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.blade_ring.extract_blade as _moved

    sys.modules[__name__] = _moved
