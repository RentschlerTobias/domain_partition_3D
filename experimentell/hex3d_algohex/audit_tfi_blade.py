"""Moved to `hexdp.quality.audit_tfi_blade` (src/hexdp/quality/audit_tfi_blade.py).

This shim keeps the old imports (`import audit_tfi_blade`) and command lines
(`python experimentell/hex3d_algohex/audit_tfi_blade.py`) working while the code lives in the hexdp package.
"""
import sys

if __name__ == "__main__":
    import runpy

    runpy.run_module("hexdp.quality.audit_tfi_blade", run_name="__main__", alter_sys=True)
else:
    import hexdp.quality.audit_tfi_blade as _moved

    sys.modules[__name__] = _moved
