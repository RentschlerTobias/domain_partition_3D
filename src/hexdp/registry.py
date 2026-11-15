"""Plugin discovery through entry points.

Every kind of plugin has one entry-point group. A package, hexdp itself or an
external one, registers a plugin in its ``pyproject.toml``::

    [project.entry-points."hexdp.strategies"]
    sweep = "hexdp.strategies.sweep:SweepStrategy"

and ``load("strategies", "sweep")`` returns an instance. Entries are loaded
lazily: listing a group imports nothing, so a plugin with a heavy runtime
(AlgoHex, dtOO, gmsh) only costs something when it is used.
"""

from dataclasses import dataclass
from importlib import import_module, metadata
from typing import Any

KINDS = (
    "detectors",
    "strategies",
    "ops",
    "discretisers",
    "fillers",
    "exporters",
    "interfaces",
)


class PluginError(LookupError):
    """A plugin is unknown or does not load."""


@dataclass(frozen=True)
class PluginEntry:
    """One registered plugin, not yet imported."""

    kind: str
    name: str
    target: str
    distribution: str | None


def group(kind: str) -> str:
    """Entry-point group name of a plugin kind."""
    if kind not in KINDS:
        raise PluginError(f"unknown plugin kind '{kind}', expected one of {', '.join(KINDS)}")
    return f"hexdp.{kind}"


def entries(kind: str) -> list[PluginEntry]:
    """All plugins registered for ``kind``, sorted by name."""
    found = []
    for ep in metadata.entry_points(group=group(kind)):
        dist = ep.dist.name if ep.dist is not None else None
        found.append(PluginEntry(kind, ep.name, ep.value, dist))
    return sorted(found, key=lambda e: e.name)


def names(kind: str) -> list[str]:
    return [e.name for e in entries(kind)]


def load(kind: str, name: str, validate: bool = True, **kwargs: Any) -> Any:
    """Import and return the plugin ``name`` of ``kind``.

    A class is instantiated with ``kwargs``; any other object (a module-level
    instance or a factory result) is returned as it is. With ``validate`` the
    object is checked against its contract (``hexdp.contracts``).
    """
    if ":" in name:                                   # direct reference "package.module:Object"
        mod, _, attr = name.partition(":")
        try:
            obj = getattr(import_module(mod), attr)
        except Exception as exc:
            raise PluginError(f"{kind} plugin '{name}' does not load: {exc}") from exc
        return _finish(kind, name, obj, validate, kwargs)
    matches = [ep for ep in metadata.entry_points(group=group(kind)) if ep.name == name]
    if not matches:
        known = ", ".join(names(kind)) or "none registered"
        raise PluginError(f"no {kind} plugin '{name}' (known: {known})")
    if len(matches) > 1:
        dists = ", ".join(sorted(str(ep.dist.name if ep.dist else "?") for ep in matches))
        raise PluginError(f"{kind} plugin '{name}' is registered more than once ({dists})")
    try:
        obj = matches[0].load()
    except Exception as exc:  # report which plugin broke, not just the import error
        raise PluginError(f"{kind} plugin '{name}' ({matches[0].value}) does not load: {exc}") from exc
    return _finish(kind, name, obj, validate, kwargs)


def _finish(kind, name, obj, validate, kwargs):
    plugin = obj(**kwargs) if isinstance(obj, type) else obj
    if validate:
        from hexdp.contracts import validate as check

        problems = check(kind, plugin)
        if problems:
            raise PluginError(f"{kind} plugin '{name}' breaks its contract: {'; '.join(problems)}")
    return plugin
