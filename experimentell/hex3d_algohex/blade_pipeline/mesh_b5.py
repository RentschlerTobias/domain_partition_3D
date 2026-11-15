"""Moved to `hexdp.strategies.blade_ring.mesh_b5` (src/hexdp/strategies/blade_ring/mesh_b5.py).

This shim keeps the old imports (`import mesh_b5`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/mesh_b5.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.blade_ring.mesh_b5", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.blade_ring.mesh_b5 as _moved

    sys.modules[__name__] = _moved
