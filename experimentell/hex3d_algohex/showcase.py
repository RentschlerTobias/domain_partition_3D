"""Moved to `hexdp.viz.showcase` (src/hexdp/viz/showcase.py).

This shim keeps the old imports (`import showcase`) and command lines
(`python experimentell/hex3d_algohex/showcase.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.viz.showcase", run_name="__main__", alter_sys=True)
else:
    import hexdp.viz.showcase as _moved

    sys.modules[__name__] = _moved
