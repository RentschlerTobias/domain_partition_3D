"""Moved to `hexdp.io.ovm` (src/hexdp/io/ovm.py).

This shim keeps the old imports (`import ovm_io`) and command lines
(`python experimentell/hex3d_algohex/ovm_io.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.io.ovm", run_name="__main__", alter_sys=True)
else:
    import hexdp.io.ovm as _moved

    sys.modules[__name__] = _moved
