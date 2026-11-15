"""Moved to `hexdp.core.block_edges` (src/hexdp/core/block_edges.py).

This shim keeps the old imports (`import block_edges`) and command lines
(`python experimentell/hex3d_algohex/block_edges.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.core.block_edges", run_name="__main__", alter_sys=True)
else:
    import hexdp.core.block_edges as _moved

    sys.modules[__name__] = _moved
