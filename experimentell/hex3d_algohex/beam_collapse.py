"""Moved to `hexdp.ops.beam_collapse` (src/hexdp/ops/beam_collapse.py).

This shim keeps the old imports (`import beam_collapse`) and command lines
(`python experimentell/hex3d_algohex/beam_collapse.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.ops.beam_collapse", run_name="__main__", alter_sys=True)
else:
    import hexdp.ops.beam_collapse as _moved

    sys.modules[__name__] = _moved
