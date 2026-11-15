"""Moved to `hexdp.quality.mesh_quality` (src/hexdp/quality/mesh_quality.py).

This shim keeps the old imports (`import mesh_quality`) and command lines
(`python experimentell/hex3d_algohex/mesh_quality.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.quality.mesh_quality", run_name="__main__", alter_sys=True)
else:
    import hexdp.quality.mesh_quality as _moved

    sys.modules[__name__] = _moved
