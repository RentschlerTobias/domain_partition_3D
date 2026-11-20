# hexdp compute backends

hexdp itself is a Python package (`pip install -e .`): gmsh, numpy, scipy, meshio, networkx and the `mcp` package
for the server. Everything else is an optional backend. `hexdp mcp` reports them with the tool `check_environment`;
a strategy whose backend is missing is listed by `list_strategies` as not applicable, with the reason.

**Rule for agents:** if a backend is missing, first check whether the available methods block the geometry (sweep,
polycube, composite building blocks such as `band_layout` / `revolve_fill` / `section_sweep`). Only if they cannot,
ask the user whether to install it - with what it would enable and what it costs (below). Never install without
consent.

| backend | found via | enables | without it |
|---|---|---|---|
| gmsh (required) | Python package | loading .msh/.step, tet meshes, `.msh` export | nothing works |
| AlgoHex | `HexMeshing` on PATH or `$HEXMESHING_BIN` | strategy `algohex`: general topology from a frame field | geometries with no sweep direction / template / polycube labeling |
| pythonOCC | `$HEXDP_OCC_PYTHON` (default `/opt/venv-dtoo/bin/python`) | strategies `sweep`, `cgrid` (exact CAD sections) | CAD sweeps |
| dtOO | `dtOOPythonSWIG` importable in the OCC python | block `dtoo_rotor` (parametric dtOO machines from `params.json`) | runner workflows from dtOO parameter files |
| evocube | `$HEXDP_EVOCUBE_BIN` (default `/opt/stack/ext/evocube/build`) | strategy `polycube` | polycube blockings (genus-0 parts) |
| mmg | `$MMG3D` | envelope cut before AlgoHex (blunt trailing edges) | AlgoHex on blunt edges may fail |
| C4HexMeshing | `$HEXDP_C4HEX_BIN` (default `/opt/stack/ext/C4HexMeshing/build/Build/bin/c4hex_cli`) | second extractor for a seamless map that HexEx rejects (`extractor: c4hex` of the algohex strategy) | AlgoHex runs whose integer-grid map is folded end without a hex mesh |

## AlgoHex - facts before using it

From the tistos datasets on this project (`/opt/stack/datasets/tistos_fluid_keepbl_2026-11-17`, 973 runs, and
`tistos_s1_2026-10-12`, 894 runs; tet target n = 2000, 8 threads, 2026-10/11):

| | |
|---|---|
| run time per geometry (successful runs) | median 35 min, 10 % under 25 min, 90 % under 76 min, longest 3 h (keep-BL) / 5.9 h (s1) |
| where the time goes (sample 0131) | frame field ~1 min, singularity optimisation ~2 min, integrability ~10 min, parametrisation the rest (up to hours) |
| no result | keep-BL 8.9 %: invalid parametrisation detected and stopped early 6.6 %, 4 h limit 2.1 %, crash 0.3 %; s1 1.6 % crashes |
| result off the surface | ~2 % of the hex meshes deviate strongly from the boundary (rejected by beam) |
| one thread | sample 0131 did not finish in 4 h |
| reproducibility | same thread count -> same hex mesh (only the entity order differs; hexdp sorts it canonically) |
| block count | varies strongly within one family (about 10 to over 200 blocks after beam); thin gaps drive it up |
| install | build from source: CoMISo, OpenVolumeMesh, HexEx, QGP3D, MC3D, TS3D, Eigen, glog, ...; typically 30-60 min with a C++17 toolchain; or the `algohex:portable` container image (see `src/hexdp/strategies/algohex/run_algohex.py`) |
| input | a tet mesh with feature tags in AlgoHex's VTK format - hexdp makes it from a `.msh` (`tet_prep`) |

There is no success guarantee: plan for about one geometry in ten without a result, and reuse an existing
`hex_hex_x.ovm` (`run_strategy(..., reuse=...)`) instead of recomputing.

## CPUs

hexdp has no CPU default of its own. It runs on the CPUs its process may use (its affinity); restrict it from
outside like any program (`taskset -c 0-7 hexdp mcp`). Over MCP every computing call (`run_strategy`,
`validate_workflow`) takes `cpus` ('8' = the first 8 available, or a list '0-7'); the agent asks the user first.
`HEXDP_MCP_CPUS` can preset the answer (unattended runs).

## Installing evocube (polycube)

GPL-3, used as external programs only. Build with the viewer off: `cmake .. -DEVOCUBE_VIEWER=OFF` then
`make prep_cli evolabel_cli polycube_withHexEx` (the local patches: `/opt/stack/ext/evocube`, see its README).
It downloads its dependencies (libigl, geogram) at configure time.

**Local patch (2026-11-20):** `app/evolabel_cli.cpp:241` and `app/evolabel.cpp:244` did
`archive.getIndiv(0)` for the final labeling, which is the *worst* individual of the archive (it is sorted by
descending score, `bestScore()` is `back()`) - changed to `archive.bestIndividual()`. Two of five `CELLS 0` models
start producing hexes with it. The build cache in `build/` still names the pre-refactor source path, so a rebuild
needs `mkdir -p /work/ext && ln -s /opt/stack/ext/evocube /work/ext/evocube` (or a fresh build dir).

**Determinism:** the label search is not reproducible with more than one thread, even though `srand(1)` is fixed:
the mutation loop is OpenMP-parallel and all threads share the global `rand()`. `OMP_NUM_THREADS=1` reproduces
(measured on n08b); a proper fix is a per-iteration RNG stream.


## C4HexMeshing - a second extractor after AlgoHex

GPL-3, used as an external program only (`src/hexdp/strategies/algohex/extract_c4hex.py`). Built here at
`/opt/stack/ext/C4HexMeshing` (`build/Build/bin/c4hex_cli`, no `LD_LIBRARY_PATH` needed; it finds its libraries by
RPATH). It is the pure extractor of Brueckler et al. 2023: motorcycle complex (MC3D) of a seamless map, optional
collapse of zero-extent cells, quantisation via QGP3D, then hex extraction and feature-preserving smoothing.

| | |
|---|---|
| input | the `.hexex` seamless map, which `run_algohex.py` writes with `--checkpoints` (`HexMeshing --sm-out-path`); the map is written before AlgoHex checks its own integer-grid map, so it exists even when that check fails |
| output | `<prefix>_<n>_hex.ovm` - OpenVolumeMesh, which `hexdp.io.ovm` and the beam read directly |
| when it helps | AlgoHex computes a seamless map but its extraction or integer-grid map is folded ("Detected N Flipped ... Tets after the Sanitization", "Generating IGM failed with error code 22"). C4Hex quantises the same map itself instead of finishing AlgoHex's IGM, so a folded IGM does not stop it |
| when it does not | a bad frame field/seamless map stays bad; and a map that even the MC-based quantisation cannot untangle still fails (exit 20-22). The `fandisk_algohex.hexex` sample of the C4Hex repo fails this way (exit 22 after 17 min) - the bundled samples are not all solvable |
| cost | 2 s for a unit cube, 17 min and 3.5 GB for the 28 MB fandisk sample; single-threaded in the parts measured so far |
| useful flags | `--collapse` (collapse the zero-arc MC - the coarse variant), `--block-structured --times-minimal K`, `--num-hexes N`, `--optimize-base-mesh`, `--untangling-iter`, `--iqp-time-limit` (seconds), `--split-selfadjacent` |

