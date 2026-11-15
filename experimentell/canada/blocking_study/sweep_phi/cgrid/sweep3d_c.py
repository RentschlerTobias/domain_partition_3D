"""Moved to `hexdp.strategies.cgrid.sweep3d_c` (src/hexdp/strategies/cgrid/sweep3d_c.py).

This shim keeps the old imports (`import sweep3d_c`) and command lines
(`python experimentell/canada/blocking_study/sweep_phi/cgrid/sweep3d_c.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.cgrid.sweep3d_c", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.cgrid.sweep3d_c as _moved

    sys.modules[__name__] = _moved
