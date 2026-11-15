"""Moved to `hexdp.strategies.algohex.tet_prep_named` (src/hexdp/strategies/algohex/tet_prep_named.py).

This shim keeps the old imports (`import tet_prep_named`) and command lines
(`python experimentell/hex3d_algohex/tet_prep_named.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.algohex.tet_prep_named", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.algohex.tet_prep_named as _moved

    sys.modules[__name__] = _moved
