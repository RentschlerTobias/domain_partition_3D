"""Moved to `hexdp.strategies.blade_ring.core_tfi` (src/hexdp/strategies/blade_ring/core_tfi.py).

This shim keeps the old imports (`import core_tfi`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/core_tfi.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.blade_ring.core_tfi", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.blade_ring.core_tfi as _moved

    sys.modules[__name__] = _moved
