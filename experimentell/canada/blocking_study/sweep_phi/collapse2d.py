"""Moved to `hexdp.strategies.sweep.collapse2d` (src/hexdp/strategies/sweep/collapse2d.py).

This shim keeps the old imports (`import collapse2d`) and command lines
(`python experimentell/canada/blocking_study/sweep_phi/collapse2d.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.sweep.collapse2d", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.sweep.collapse2d as _moved

    sys.modules[__name__] = _moved
