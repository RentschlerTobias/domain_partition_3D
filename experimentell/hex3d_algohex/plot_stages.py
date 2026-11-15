"""Moved to `hexdp.viz.plot_stages` (src/hexdp/viz/plot_stages.py).

This shim keeps the old imports (`import plot_stages`) and command lines
(`python experimentell/hex3d_algohex/plot_stages.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.viz.plot_stages", run_name="__main__", alter_sys=True)
else:
    import hexdp.viz.plot_stages as _moved

    sys.modules[__name__] = _moved
