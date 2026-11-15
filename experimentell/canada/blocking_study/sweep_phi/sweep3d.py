"""Moved to `hexdp.strategies.sweep.sweep3d` (src/hexdp/strategies/sweep/sweep3d.py).

This shim keeps the old imports (`import sweep3d`) and command lines
(`python experimentell/canada/blocking_study/sweep_phi/sweep3d.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.sweep.sweep3d", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.sweep.sweep3d as _moved

    sys.modules[__name__] = _moved
