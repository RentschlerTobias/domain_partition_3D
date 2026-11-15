"""Moved to `hexdp.quality.check_mesh_openfoam` (src/hexdp/quality/check_mesh_openfoam.py).

This shim keeps the old imports (`import check_mesh_openfoam`) and command lines
(`python experimentell/hex3d_algohex/datagen/check_mesh_openfoam.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.quality.check_mesh_openfoam", run_name="__main__", alter_sys=True)
else:
    import hexdp.quality.check_mesh_openfoam as _moved

    sys.modules[__name__] = _moved
