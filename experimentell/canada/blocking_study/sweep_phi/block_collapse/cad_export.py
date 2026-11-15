"""Moved to `hexdp.io.cad_export` (src/hexdp/io/cad_export.py).

This shim keeps the old imports (`import cad_export`) and command lines
(`python experimentell/canada/blocking_study/sweep_phi/block_collapse/cad_export.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.io.cad_export", run_name="__main__", alter_sys=True)
else:
    import hexdp.io.cad_export as _moved

    sys.modules[__name__] = _moved
