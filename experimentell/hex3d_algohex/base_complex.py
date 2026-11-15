"""Moved to `hexdp.core.base_complex` (src/hexdp/core/base_complex.py).

This shim keeps the old imports (`import base_complex`) and command lines
(`python experimentell/hex3d_algohex/base_complex.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.core.base_complex", run_name="__main__", alter_sys=True)
else:
    import hexdp.core.base_complex as _moved

    sys.modules[__name__] = _moved
