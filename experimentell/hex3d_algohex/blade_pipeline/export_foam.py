"""Moved to `hexdp.export.foam` (src/hexdp/export/foam.py).

This shim keeps the old imports (`import export_foam`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/export_foam.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.export.foam", run_name="__main__", alter_sys=True)
else:
    import hexdp.export.foam as _moved

    sys.modules[__name__] = _moved
