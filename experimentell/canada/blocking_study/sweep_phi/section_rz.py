"""Moved to `hexdp.strategies.sweep.section_rz` (src/hexdp/strategies/sweep/section_rz.py).

This shim keeps the old imports (`import section_rz`) and command lines
(`python experimentell/canada/blocking_study/sweep_phi/section_rz.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.sweep.section_rz", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.sweep.section_rz as _moved

    sys.modules[__name__] = _moved
