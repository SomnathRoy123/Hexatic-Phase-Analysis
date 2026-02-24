# Translational Order (Standalone)

This directory is a standalone translational-order analysis pipeline split from the hexatic-only workflow.

## What it computes

For each frame (`time_<idx>.dat`):
1. Cluster particles using bond cutoff `lbond`.
2. Compute cluster COMs.
3. Build Delaunay neighbors on COMs.
4. Compute local `psi6` and frame-wise global orientation angle `theta_G`.
5. Accumulate:
   - Radial distribution function `g(r)` (`gr_avg_time_<start>_<end>.dat`)
   - Translational correlation `gT(r)` (`gt_avg_time_<start>_<end>.dat`)

## Build

```bash
make
```

## Run

```bash
./translational_order DATA_DIR START_IDX END_IDX OUT_DIR LBOND DR A_LATTICE USE_PBC BOX_X BOX_Y
```

Example:

```bash
./translational_order ./data/ 1000 2000 ./out/ 1.5 0.5 1.12 1 180.0 180.0
```
