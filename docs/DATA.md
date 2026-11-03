# Data inventory (block-generation work)

Status: 2026-11-03. Everything below lives outside the git repository; this file is the
manifest that makes it findable again.

## 1. tistos_s1 corpus (the ~1000 designs)

| Field | Value |
|---|---|
| Archive | `tistos_s1_2026-10-12.tar.xz` (was `datagen_tistos_s1_papaya_2026-10-12.tar.xz`) |
| Location | `/opt/stack/datasets/tistos_s1_2026-10-12/` |
| Size | 16,703,213,184 bytes (16 GiB) |
| SHA256 | `d83405c4525e0b66ad1ecf7a0a8fb1e9565f86be65ff6883d4f8223ba4af4273` |
| Top-level in the archive | `datagen/` (runs/, locks/, lane markers, pack.sh) |

Verify and extract:

    cd /opt/stack/datasets/tistos_s1_2026-10-12
    sha256sum -c tistos_s1_2026-10-12.tar.xz.sha256
    tar -xf tistos_s1_2026-10-12.tar.xz -C "$DP3D_DATA_ROOT"

Contents: 1000 run directories `tistos_s1_0000 … tistos_s1_0999`.

| Quantity | Count | Meaning |
|---|---|---|
| run directories | 1000 | drawn Sobol designs |
| `hex_hex_x*.ovm` | 880 | AlgoHex block structure present |
| `sample.npz` | 862 | finished sample (12 blocks / 12 cuboids, meshtron keys) |

A finished run contains `params.json`, `case/mesh.msh` (exact dtOO geometry),
`tet.vtk`, `hex_hex_x_16c.ovm`, `sample.npz` and the stage logs.

## 2. keep-BL case (not in the archive)

`/work/datagen_keepbl/runs/tistos_s1_0131` (388 MB) — the tistos case with the blade
ring, the TFI, the audits and the visualisations produced on 2026-11-02/03. Created
after the archive was packed, so it must be copied separately (it is small).

Persistent copy: `/opt/stack/datasets/keepbl/tistos_s1_0131/` (387 MB, 2026-11-03).

## 3. Other data sets

| Set | Location | Size | Notes |
|---|---|---|---|
| `blade_blocking_20261021` | `/opt/stack/datasets/blade_blocking_20261021/` | 30 MB | blade-blocking deliverable (`keep_bl/`, `blocking_0131_gap_opt/`) |
| `meshtron_flow_20261012` | `/opt/stack/datasets/meshtron_flow_20261012/` | 3.7 GB | transport archive of the flow-matching work; already extracted into `meshtron/data` + `meshtron/experiments/flow_topology` |
| `courier_datagen_2026-10-06` | `/opt/stack/archive/courier_datagen_2026-10-06/` | 1.3 GB | aborted 400-geometry job, 11 finished samples (all included in the corpus above) |

## 3b. Case archives (inputs, runs, results)

Copied from the sandbox on 2026-11-03; each folder has its own `README.md` with the
layout and a "what is input / what is derived" section.

| Folder | Size | Contents |
|---|---|---|
| `/opt/stack/cases/canada/` | 24 GB | canada case `E1_12685`: `cad/` and `case/` (inputs), `dt/` variants, `pilot/` and `sweep/` runs with `metrics_*.json`, `envelope/` experiments, `vtk/` renders |
| `/opt/stack/cases/calib/` | 7.5 GB | calibration: `report_hexme.json`, `report_charite*.json`, `hexme/meshes/` (190 VTK meshes), `sweep/`, `hexme_sweep/`, `smoke/`, `collapse_probe/` |

## 4. Working copies (ephemeral, container-local /work)

| Path | Size | Reproducible from |
|---|---|---|
| `/work/datagen/runs/` | 103 GB | the archive in section 1 |
| `/work/datagen_keepbl/runs/tistos_s1_0131/` | 388 MB | section 2 (this is the only original) |
| `/work/canada/` | 24 GB | canada pipeline (`gen_mesh_canada.py` etc.) |
| `/work/calib/` | 7.5 GB | calibration sweep, `calib_sweep.py` |
| `/work/tools/` | 738 MB | mmg + algohex-size builds |

`/work` is **container-local and disappears with the container**. Everything in sections
1-3b has a persistent copy under `/opt/stack`; the `/work` entries are just the working
copies. The only exception is `/work/tools/` (tool builds, not data).

## 5. Visualisations

Curated, persistent: `/opt/stack/viz/` (scene + PNG + provenance JSON per item, with an
`index.md`). Generator output stays in `/work/...` and is not curated.
