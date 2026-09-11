# Decision log — dataset sampler strategy (grilling session 2026-09-11)

Topic: shape the hex3d dataset generation so it is optimizer-ready in
structure, without wiring the optimizer now.

## Live design tree

```
Generate machine geometry dataset (optimizer-shaped, optimizer NOT wired)
├── [x] Sampling scheme: Sobol with --strategy seam      (Q1)
├── [ ] Sampler state/persistence contract               (next)
├── [ ] First batch size (blocked on T10 throughput + yield)
├── [x] n per machine = 2 (n=2000, n=8000)               (T6, settled)
├── [x] Runner: batch_samples.slurm with CANDS           (exists)
├── [x] geometry export: DtooAdapter.export_mesh() in dtOO enroot container (T14)
└── [x] no optimizer connection in this phase            (user constraint)
```

## Q1 — Sampling scheme

**Decision:** Sobol (scipy.stats.qmc.Sobol) with a `--strategy` flag on the
sampler script, so a future Differential-Evolution strategy sweeps in as an
added plugin that fills the same candidate layout.

**Options considered:**
- (a) Sobol
- (b) Latin Hypercube
- (c) uniform random
- (d) Sobol + explicit `--strategy` seam (CHOSEN — identical to (a) but names
  the seam the optimizer will use later)

**Why:** In 30 dimensions uniform draws clump; Sobol fillets the box evenly at
any N, and is incremental: with a fixed seed, drawing additional points
continues the sequence and fills the gaps left by an earlier job — a job
killed at a walltime cap can simply be resubmitted and the resulting set
behaves like one planned sample. LHS cannot be extended without breaking its
distribution property; uniform random leaves real holes at 30-D. The later
DE strategy starts its population behind the last dataset index, so the
sampler is the single seam between "data generation now" and
"optimization-relevant".

## Practical constraints recorded

- Sobol balance is best at 2^m sample counts: round the target up to the next
  power of two and allow job boundaries between chunks.
- Resumability needs persistent sampler state (which draw indices are already
  materialized), otherwise a resubmit replays from scratch. To be settled in
  the next question.
- Unknown dtOO yield rate (DTOO_FAIL_PENALTY exists) multiplies onto the
  target count; measure on the first batch, do not assume (already recorded
  in T14 of DATASET_PIPELINE.md).

## Q2 — Sampler state for resumability

**Decision:** Option (c). A `sampler.json` metadata file beside the output
root (seed, strategy, design-bounds hash, target count) plus directory
presence as the idempotency marker. Failed draws go to the back of the queue
and are retried on the next run — the same pattern as `sample_one.sh`
("resubmit and the finished ones are skipped"), lifted to the geometry layer.

**Options considered:** (a) presence-only, (b) state file only, (c) both.

**Why:** (a) alone breaks the Sobol sequence once a draw fails mid-forest
(indices shift, duplicates appear); (b) alone creates a second source of
truth that can desync from the filesystem. (c) keeps metadata (seed, bounds
hash, target count) and progress in separate channels, neither of which can
silently contradict the data. Seed + bounds hash make later extension exact.

## Q3 — Dataset namespace

**Decision:** Option (b). New dataset root `data/dataset/sobol/machine_00NN/`
with `sampler.json` beside it. The existing 21 candidates stay under
`data/random_tistos/investigated/` as a reference pool.

**Why:** The 21 were produced by the earlier random_tistos generation AND
mostly hold tet meshes only — most have no hex-block samples yet. Mixing a
seeded Sobol stream (with seed, bounds hash, target count) into a filtered
reference pool would blur the meaning of the candidate index. In the new
root, `machine_00NN` = draw NN of the Sobol sequence is load-bearing
 addressing that the optimizer can use later; the desired directory
layout is what the optimizer phase sees as one platform.

**Correction recorded:** T6's "27 min both" dataset pair (n=2000/8000) was
settled on the T1_9 fixture, not on the 21 candidates. Not all of the 21
candidates have hex meshes; they are primarily tet-mesh suppliers.

## Q4 — Sampler code ownership

**Decision:** Option (a). The sampler driver lives in this repo
(`experimentell/hex3d_algohex/scripts/generate_machines.py`, owns layout,
sampler.json, --strategy) and reaches the dtOO world through ONE import seam:
`DtooAdapter` and `run_dtoo_export` from the sibling `eigenfrequencies` repo
via a path insert, mounted into the container as documented in
`cluster/enroot_dtoo_import.md`.

**User's own words (decisive reason):** the domain-partition pipeline will
later be integrated INTO the eigenfrequencies optimization pipeline; since
that is the destination, we already use that repo's container/import today
and do not copy the adapter into a second place.

**Rejected:** sampler in eigenfrequencies (dataset-layout knowledge in the
wrong repo); copying the adapter (duplicates documented P0 container traps:
double env.sh sourcing, no set -u, container name-mismatch hazard).
