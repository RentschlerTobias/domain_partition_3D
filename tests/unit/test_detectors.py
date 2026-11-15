import math

import numpy as np
import pytest

from hexdp.contracts import resolve_params, validate
from hexdp.detect.rotation import RotationAxisDetector, SweepDetector, fit_axis, line_distances


@pytest.fixture(scope="module")
def shapes(tmp_path_factory):
    """quarter annulus about z, the same with a box cut out (a 'blade'), a cylinder on a tilted axis, a box"""
    gmsh = pytest.importorskip("gmsh")
    d = tmp_path_factory.mktemp("shapes")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    def sector():
        s = gmsh.model.occ.addRectangle(0.5, 0, 0, 1.0, 2.0)
        gmsh.model.occ.rotate([(2, s)], 0, 0, 0, 1, 0, 0, math.pi / 2)
        gmsh.model.occ.revolve([(2, s)], 0, 0, 0, 0, 0, 1, math.pi / 2)

    def blade():
        sector()
        b = gmsh.model.occ.addBox(0.9, 0.2, 0.5, 0.4, 0.05, 0.8)
        gmsh.model.occ.cut([(3, 1)], [(3, b)])

    builds = {"sector": sector, "blade": blade,
              "tilted": lambda: gmsh.model.occ.addCylinder(1, 0, 0, 1, 1, 0, 0.4),
              "box": lambda: gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)}
    for name, build in builds.items():
        gmsh.clear()
        build()
        gmsh.model.occ.synchronize()
        gmsh.write(str(d / f"{name}.step"))
    gmsh.finalize()
    from hexdp.io.gmsh_model import load
    return {n: load(d / f"{n}.step") for n in builds}


def run(det, geo):
    return det.detect(geo, resolve_params(det.params, {}), None)


def test_contracts():
    assert validate("detectors", RotationAxisDetector()) == [] and validate("detectors", SweepDetector()) == []


def test_axis_of_a_line_fit():
    """points on a cylinder r=1 about the line through (2, 0, 0) along y, with their normals"""
    t = np.linspace(0, 2 * np.pi, 40, endpoint=False)
    y = np.linspace(-1, 1, 7)
    T, Y = np.meshgrid(t, y)
    N = np.stack([np.cos(T), np.zeros_like(T), np.sin(T)], -1).reshape(-1, 3)
    X = np.array([2.0, 0, 0]) + np.stack([np.cos(T), Y, np.sin(T)], -1).reshape(-1, 3)
    p, a, d = fit_axis(X, N, np.ones(len(X)))
    assert abs(abs(a[1]) - 1) < 1e-9 and np.allclose(p, [2, 0, 0], atol=1e-9) and d.max() < 1e-9
    assert line_distances(X, N, np.array([2.5, 0, 0]), np.array([0, 1.0, 0])).max() > 0.4


def test_sector_axis_and_sweep(shapes):
    ax = run(RotationAxisDetector(), shapes["sector"])[0]
    assert np.allclose(np.abs(ax.data["direction"]), [0, 0, 1], atol=1e-3)
    assert np.linalg.norm(ax.data["point"]) < 1e-2 and abs(ax.data["angular_span_deg"] - 90) < 0.5
    assert len(ax.data["sector_planes"]) == 2 and not ax.data["other_surfaces"]
    sw = run(SweepDetector(), shapes["sector"])
    assert [f.kind for f in sw] == ["sweep"] and sw[0].confidence > 0.9


def test_blade_breaks_the_sweep(shapes):
    ax = run(RotationAxisDetector(), shapes["blade"])[0]
    assert ax.data["other_surfaces"]                      # the blade faces are not rotational
    sw = run(SweepDetector(), shapes["blade"])
    assert [f.kind for f in sw] == ["sweep_rejected"] and "change with phi" in " ".join(sw[0].reasons)


def test_tilted_axis(shapes):
    ax = run(RotationAxisDetector(), shapes["tilted"])[0]
    assert np.allclose(ax.data["direction"], [2 ** -0.5, 2 ** -0.5, 0], atol=1e-3)
    assert np.allclose(ax.data["point"], [0.5, -0.5, 0], atol=1e-2)


def test_box_has_no_axis(shapes):
    """two parallel planes do not fix an axis: no rotation axis, no sweep"""
    assert run(RotationAxisDetector(), shapes["box"]) == [] and run(SweepDetector(), shapes["box"]) == []
