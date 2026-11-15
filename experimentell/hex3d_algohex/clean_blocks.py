"""Moved to `hexdp.core.clean_blocks` (src/hexdp/core/clean_blocks.py).

This shim keeps the old imports (`import clean_blocks`) and command lines
(`python experimentell/hex3d_algohex/clean_blocks.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.core.clean_blocks", run_name="__main__", alter_sys=True)
else:
    import hexdp.core.clean_blocks as _moved

    sys.modules[__name__] = _moved
