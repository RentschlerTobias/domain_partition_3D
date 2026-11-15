"""Moved to `hexdp.strategies.cgrid.wall_faces` (src/hexdp/strategies/cgrid/wall_faces.py).

This shim keeps the old imports (`import wall_faces`) and command lines
(`python experimentell/canada/blocking_study/sweep_phi/cgrid/wall_faces.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.strategies.cgrid.wall_faces", run_name="__main__", alter_sys=True)
else:
    import hexdp.strategies.cgrid.wall_faces as _moved

    sys.modules[__name__] = _moved
