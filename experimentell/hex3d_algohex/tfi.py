"""Moved to `hexdp.fill.tfi` (src/hexdp/fill/tfi.py).

This shim keeps the old imports (`import tfi`) and command lines
(`python experimentell/hex3d_algohex/tfi.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.fill.tfi", run_name="__main__", alter_sys=True)
else:
    import hexdp.fill.tfi as _moved

    sys.modules[__name__] = _moved
