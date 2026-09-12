# Cluster job logs

Cluster jobs write one log per job into this directory and the logs are
committed to git. This mirrors the convention in the sibling eigenfrequencies
repo (`cluster/logs/<jobid>/`), so a failed run leaves its evidence somewhere
greppable instead of only in a scratch directory that gets rotated away.

- `batch_dtoo_export_<jobid>.log` — T14 pass 2, dtOO `export_mesh` per machine.
- `batch_samples_<jobid>.log` — T10 hexablock batch (the driver tees its output
  here; the SLURM `.out` in the submit directory holds the same text).

The SLURM `.out` in the submit directory still exists, but it is not durable;
this directory is. After a run:

```bash
git add experimentell/hex3d_algohex/logs
git commit -m "cluster: dtOO export job <jobid> log"
```

What to read out of a log, in order: the guard lines (a missing image or
checkout aborts in seconds and the log says which), the container line
(`unpacking` once ever, then `already unpacked`), the yield summary at the end
(exported / failed / not attempted), and for a `batch_samples` log the
throughput block. The yield number and the throughput number together settle
the dataset size, the last open decision in DATASET_PIPELINE.md.
