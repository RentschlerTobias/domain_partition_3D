"""Moved to `hexdp.quality.nonortho` (src/hexdp/quality/nonortho.py).

This shim keeps the old imports (`import nonortho`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/nonortho.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.quality.nonortho", run_name="__main__", alter_sys=True)
else:
    import hexdp.quality.nonortho as _moved

    sys.modules[__name__] = _moved
