import importlib

import pytest

SUBPACKAGES = ["io", "core", "detect", "strategies", "ops", "discretise", "fill",
               "quality", "export", "interfaces", "pipeline", "mcp"]


@pytest.mark.parametrize("name", SUBPACKAGES)
def test_subpackage_imports_and_says_what_it_is(name):
    mod = importlib.import_module(f"hexdp.{name}")
    assert mod.__doc__ and len(mod.__doc__) > 20
