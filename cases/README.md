# Cases

One YAML per blocking job; run with `hexdp run <case>.yaml`, several with `hexdp eval cases/*.yaml`.
Input data are not in the repository: the cases read them from `$HEXDP_DATA`:

    $HEXDP_DATA/tistos_s1_0131/case/mesh.msh      dtOO hybrid mesh of tistos sample 0131
    $HEXDP_DATA/tistos_s1_0131/hex_hex_x.ovm      its AlgoHex hex mesh (reused: AlgoHex takes ~40 min)
    $HEXDP_DATA/canada/domain.brep                canada draft-tube domain solid (baffle domain)
    $HEXDP_DATA/canada/cad_blk.stp                canada draft-tube CAD with the thin walls

Workspaces go to `runs/<name>` next to the case file (ignored by git). The `reference` numbers are the
results of the legacy scripts before the migration (bit-identical regression, see tests/regression).
