"""Moved to `hexdp.quality.hexq` (src/hexdp/quality/hexq.py).

This shim keeps the old imports (`import quality`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/quality.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.quality.hexq", run_name="__main__", alter_sys=True)
else:
    import hexdp.quality.hexq as _moved

    sys.modules[__name__] = _moved
