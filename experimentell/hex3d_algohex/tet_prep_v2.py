"""Moved to `hexdp.strategies.algohex.tet_prep_v2` (src/hexdp/strategies/algohex/tet_prep_v2.py).

This shim keeps the old imports (`import tet_prep_v2`) and command lines
(`python experimentell/hex3d_algohex/tet_prep_v2.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.algohex.tet_prep_v2", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.algohex.tet_prep_v2 as _moved

    sys.modules[__name__] = _moved
