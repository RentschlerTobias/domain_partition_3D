"""Moved to `hexdp.ops.sheet_beam` (src/hexdp/ops/sheet_beam.py).

This shim keeps the old imports (`import sheet_beam`) and command lines
(`python experimentell/hex3d_algohex/sheet_beam.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.ops.sheet_beam", run_name="__main__", alter_sys=True)
else:
    import hexdp.ops.sheet_beam as _moved

    sys.modules[__name__] = _moved
