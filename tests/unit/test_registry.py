from importlib import metadata

import pytest

from hexdp import registry


class _Dist:
    name = "fake-dist"


def _ep(name, value):
    ep = metadata.EntryPoint(name=name, value=value, group="hexdp.strategies")
    return ep._for(_Dist()) if hasattr(ep, "_for") else ep


@pytest.fixture
def fake_entry_points(monkeypatch):
    """Replace the entry-point table with a few controlled entries."""
    table = {
        "hexdp.strategies": [
            _ep("dummy", "tests.unit.plugins_fixture:DummyStrategy"),
            _ep("broken", "tests.unit.plugins_fixture:DoesNotExist"),
            _ep("twice", "tests.unit.plugins_fixture:DummyStrategy"),
            _ep("twice", "tests.unit.plugins_fixture:DummyStrategy"),
            _ep("sloppy", "tests.unit.plugins_fixture:NoParams"),
        ]
    }
    monkeypatch.setattr(registry.metadata, "entry_points", lambda group: table.get(group, []))


def test_every_kind_has_a_group():
    assert [registry.group(k) for k in registry.KINDS] == [f"hexdp.{k}" for k in registry.KINDS]


def test_unknown_kind_is_rejected():
    with pytest.raises(registry.PluginError, match="unknown plugin kind"):
        registry.group("meshers")


def test_entries_are_listed_without_importing(fake_entry_points):
    assert registry.names("strategies") == ["broken", "dummy", "sloppy", "twice", "twice"]
    assert registry.names("ops") == []


def test_load_instantiates_a_class_with_kwargs(fake_entry_points):
    plugin = registry.load("strategies", "dummy", tol=0.1)
    assert plugin.name == "dummy" and plugin.tol == 0.1


def test_unknown_name_lists_the_known_ones(fake_entry_points):
    with pytest.raises(registry.PluginError, match="known: broken, dummy, sloppy, twice"):
        registry.load("strategies", "sweep")


def test_broken_plugin_names_itself(fake_entry_points):
    with pytest.raises(registry.PluginError, match="'broken'.*does not load"):
        registry.load("strategies", "broken")


def test_duplicate_registration_is_an_error(fake_entry_points):
    with pytest.raises(registry.PluginError, match="more than once"):
        registry.load("strategies", "twice")


def test_contract_is_checked_on_load(fake_entry_points):
    with pytest.raises(registry.PluginError, match="breaks its contract.*params.*'run'"):
        registry.load("strategies", "sloppy")
    assert registry.load("strategies", "sloppy", validate=False).name == "noparams"
