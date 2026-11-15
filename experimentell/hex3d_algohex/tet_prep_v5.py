"""Moved to `hexdp.strategies.algohex.tet_prep` (src/hexdp/strategies/algohex/tet_prep.py).

This shim keeps the old imports (`import tet_prep_v5`) and command lines
(`python experimentell/hex3d_algohex/tet_prep_v5.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.algohex.tet_prep", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.algohex.tet_prep as _moved

    sys.modules[__name__] = _moved
