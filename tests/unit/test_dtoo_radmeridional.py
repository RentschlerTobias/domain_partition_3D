"""dtoo_radmeridional: the builder import survives a dtOO built without some SWIG classes (no dtOO needed)."""
import sys
import types

import pytest

from hexdp.blocks import dtoo_radmeridional as dr


def test_import_builders_adds_placeholders(monkeypatch, tmp_path):
    swig = types.ModuleType("dtOOPythonSWIG")
    monkeypatch.setitem(sys.modules, "dtOOPythonSWIG", swig)
    pkg = tmp_path / "dtOOPythonApp" / "builder"
    pkg.mkdir(parents=True)
    (tmp_path / "dtOOPythonApp" / "__init__.py").write_text("")
    (pkg / "__init__.py").write_text("from dtOOPythonSWIG import ofOpenFOAMCase\nfrom dtOOPythonSWIG import ofOther\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    for k in [k for k in sys.modules if k.startswith("dtOOPythonApp")]:
        monkeypatch.delitem(sys.modules, k)
    assert dr.import_builders() == ["ofOpenFOAMCase", "ofOther"]
    with pytest.raises(RuntimeError, match="not part of this dtOO build"):
        swig.ofOpenFOAMCase()


def test_plain():
    assert dr.plain({"a": (1, 2.0), "b": object}) == {"a": [1, 2.0], "b": str(object)}
