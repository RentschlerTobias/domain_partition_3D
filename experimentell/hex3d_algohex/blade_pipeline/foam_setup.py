"""Moved to `hexdp.export.foam_setup` (src/hexdp/export/foam_setup.py).

This shim keeps the old imports (`import foam_setup`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/foam_setup.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.export.foam_setup", run_name="__main__", alter_sys=True)
else:
    import hexdp.export.foam_setup as _moved

    sys.modules[__name__] = _moved
