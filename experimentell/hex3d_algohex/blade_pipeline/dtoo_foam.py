"""Moved to `hexdp.export.dtoo_foam` (src/hexdp/export/dtoo_foam.py).

This shim keeps the old imports (`import dtoo_foam`) and command lines
(`python experimentell/hex3d_algohex/blade_pipeline/dtoo_foam.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.export.dtoo_foam", run_name="__main__", alter_sys=True)
else:
    import hexdp.export.dtoo_foam as _moved

    sys.modules[__name__] = _moved
