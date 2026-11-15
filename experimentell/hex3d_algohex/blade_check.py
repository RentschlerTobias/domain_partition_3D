"""Moved to `hexdp.quality.blade_check` (src/hexdp/quality/blade_check.py).

This shim keeps the old imports (`import blade_check`) and command lines
(`python experimentell/hex3d_algohex/blade_check.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.quality.blade_check", run_name="__main__", alter_sys=True)
else:
    import hexdp.quality.blade_check as _moved

    sys.modules[__name__] = _moved
