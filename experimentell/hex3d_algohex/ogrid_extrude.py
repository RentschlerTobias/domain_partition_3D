"""Moved to `hexdp.ops.ogrid_extrude` (src/hexdp/ops/ogrid_extrude.py).

This shim keeps the old imports (`import ogrid_extrude`) and command lines
(`python experimentell/hex3d_algohex/ogrid_extrude.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.ops.ogrid_extrude", run_name="__main__", alter_sys=True)
else:
    import hexdp.ops.ogrid_extrude as _moved

    sys.modules[__name__] = _moved
