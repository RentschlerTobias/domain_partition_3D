"""Moved to `hexdp.export.foam_compare` (src/hexdp/export/foam_compare.py).

This shim keeps the old imports (`import foam_compare`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/foam_compare.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.export.foam_compare", run_name="__main__", alter_sys=True)
else:
    import hexdp.export.foam_compare as _moved

    sys.modules[__name__] = _moved
