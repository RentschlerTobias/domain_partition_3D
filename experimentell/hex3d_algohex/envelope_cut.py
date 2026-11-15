"""Moved to `hexdp.strategies.algohex.envelope_cut` (src/hexdp/strategies/algohex/envelope_cut.py).

This shim keeps the old imports (`import envelope_cut`) and command lines
(`python experimentell/hex3d_algohex/envelope_cut.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.algohex.envelope_cut", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.algohex.envelope_cut as _moved

    sys.modules[__name__] = _moved
