"""Moved to `hexdp.strategies.algohex.run_algohex` (src/hexdp/strategies/algohex/run_algohex.py).

This shim keeps the old imports (`import run_algohex`) and command lines
(`python experimentell/hex3d_algohex/run_algohex.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.algohex.run_algohex", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.algohex.run_algohex as _moved

    sys.modules[__name__] = _moved
