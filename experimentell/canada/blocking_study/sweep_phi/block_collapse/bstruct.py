"""Moved to `hexdp.ops.bstruct` (src/hexdp/ops/bstruct.py).

This shim keeps the old imports (`import bstruct`) and command lines
(`python experimentell/canada/blocking_study/sweep_phi/block_collapse/bstruct.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.ops.bstruct", run_name="__main__", alter_sys=True)
else:
    import hexdp.ops.bstruct as _moved

    sys.modules[__name__] = _moved
