"""Moved to `hexdp.core.block_faces` (src/hexdp/core/block_faces.py).

This shim keeps the old imports (`import block_faces`) and command lines
(`python experimentell/hex3d_algohex/block_faces.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.core.block_faces", run_name="__main__", alter_sys=True)
else:
    import hexdp.core.block_faces as _moved

    sys.modules[__name__] = _moved
