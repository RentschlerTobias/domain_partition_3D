"""Moved to `hexdp.fill.curve_postprocess` (src/hexdp/fill/curve_postprocess.py).

This shim keeps the old imports (`import curve_postprocess`) and command lines
(`python experimentell/hex3d_algohex/curve_postprocess.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.fill.curve_postprocess", run_name="__main__", alter_sys=True)
else:
    import hexdp.fill.curve_postprocess as _moved

    sys.modules[__name__] = _moved
