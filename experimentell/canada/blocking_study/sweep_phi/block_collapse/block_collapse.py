"""Moved to `hexdp.ops.block_collapse` (src/hexdp/ops/block_collapse.py).

This shim keeps the old imports (`import block_collapse`) and command lines
(`python experimentell/canada/blocking_study/sweep_phi/block_collapse/block_collapse.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.ops.block_collapse", run_name="__main__", alter_sys=True)
else:
    import hexdp.ops.block_collapse as _moved

    sys.modules[__name__] = _moved
