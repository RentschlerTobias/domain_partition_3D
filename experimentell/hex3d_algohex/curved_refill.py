"""Moved to `hexdp.fill.curved_refill` (src/hexdp/fill/curved_refill.py).

This shim keeps the old imports (`import curved_refill`) and command lines
(`python experimentell/hex3d_algohex/curved_refill.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.fill.curved_refill", run_name="__main__", alter_sys=True)
else:
    import hexdp.fill.curved_refill as _moved

    sys.modules[__name__] = _moved
