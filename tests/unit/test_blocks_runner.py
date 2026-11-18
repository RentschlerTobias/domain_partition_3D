"""Building blocks band_layout + revolve_fill + section_sweep composed into a solid rotor, on a synthetic rotor (no
dtOO): hub cylinder + 4 staggered elliptic blades with a fillet that widens the footprint on the hub."""
import numpy as np
import pytest

from hexdp.blocks import compose
from hexdp.blocks.band import band_layout
from hexdp.blocks.revolve import revolve_fill
from hexdp.blocks.rotor import rotor_description
from hexdp.blocks.sweep import section_sweep


def synthetic_rotor(nb=4, r_hub=0.5, r_fil=0.6, r_tip=1.4, nu=41, nw=9):
    t_cut = np.linspace(0, 2 * np.pi, 7)                 # 6 patches around the profile
    st = np.radians(55)                                  # stagger
    xc, zc = 0.0, 1.25

    def ring(r, a, b):
        out = []
        for k in range(6):
            t = np.linspace(t_cut[k], t_cut[k + 1], nu)
            x0, z0 = a * np.cos(t), b * np.sin(t)
            x = xc + x0 * np.cos(st) - z0 * np.sin(st)
            z = zc + x0 * np.sin(st) + z0 * np.cos(st)
            phi = x / r_hub
            out.append(np.c_[r * np.cos(phi), r * np.sin(phi), z])
        return out

    def patches(r0, r1, a0, a1, b0, b1):
        W = np.linspace(0, 1, nw)
        rings = [ring(r0 + w * (r1 - r0), a0 + w * (a1 - a0), b0 + w * (b1 - b0)) for w in W]
        return np.array([np.stack([rings[j][k] for j in range(nw)], axis=1) for k in range(6)])
    return dict(r_hub=r_hub, z0=0.0, z1=2.5, nb=nb,
                fillet=patches(r_hub, r_fil, 0.36, 0.3, 0.14, 0.08), blade=patches(r_fil, r_tip, 0.3, 0.3, 0.08, 0.08))


@pytest.mark.parametrize("fillet, r_bore", [(True, 0.0), (True, 0.15), (False, 0.0), (False, 0.15)])
def test_synthetic_rotor_is_valid_and_closed(fillet, r_bore):
    rot = rotor_description(synthetic_rotor(), fillet=fillet)
    band = band_layout(rot, h=0.06, core_pairing=r_bore == 0)
    hub = revolve_fill(band, r_bore=r_bore)
    blades = section_sweep(band, rot["stages"], h=0.06)
    bs, info = compose([hub, blades])
    s = bs.summary()
    assert s["valid"], s["issues"]
    assert bs.boundary_topology() == {"chi": 2 if r_bore == 0 else 0, "components": 1}   # a bore: a torus
    names = {bs.meta["surface_names"][fi.surface] for fi in bs.face_info.values()}
    assert names == {"hub_z0", "hub_z1", "hub_mantle", "blade", "blade_tip"} | ({"hub_bore"} if r_bore else set())
