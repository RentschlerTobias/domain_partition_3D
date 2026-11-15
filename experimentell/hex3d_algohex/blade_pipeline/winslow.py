"""Moved to `hexdp.fill.winslow` (src/hexdp/fill/winslow.py).

This shim keeps the old imports (`import winslow`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/winslow.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.fill.winslow", run_name="__main__", alter_sys=True)
else:
    import hexdp.fill.winslow as _moved

    sys.modules[__name__] = _moved
