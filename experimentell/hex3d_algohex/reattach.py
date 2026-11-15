"""Moved to `hexdp.fill.reattach` (src/hexdp/fill/reattach.py).

This shim keeps the old imports (`import reattach`) and command lines
(`python experimentell/hex3d_algohex/reattach.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.fill.reattach", run_name="__main__", alter_sys=True)
else:
    import hexdp.fill.reattach as _moved

    sys.modules[__name__] = _moved
