"""Moved to `hexdp.io.vtk_export` (src/hexdp/io/vtk_export.py).

This shim keeps the old imports (`import export_vtk`) and command lines
(`python experimentell/hex3d_algohex/export_vtk.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.io.vtk_export", run_name="__main__", alter_sys=True)
else:
    import hexdp.io.vtk_export as _moved

    sys.modules[__name__] = _moved
