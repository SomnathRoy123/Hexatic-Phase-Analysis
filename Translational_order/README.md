# Translational Order (Standalone)

This directory is a standalone translational-order analysis pipeline split from the hexatic-only workflow.

## What it computes

For each frame (`time_<idx>.dat`) independently (no cluster identity tracking across frames):
1. Cluster particles using bond cutoff `lbond`.
2. Compute cluster COMs for that frame.
3. Build Delaunay neighbors on COMs.
4. Compute local `psi6` and frame-wise global orientation angle `theta_G`.
5. Accumulate outputs:
   - Radial distribution function `g(r)` from COM pairs (`gr_avg_time_<start>_<end>.dat`).
   - Translational correlation `gT(r)` (`gt_avg_time_<start>_<end>.dat`).

## Lattice constant for reciprocal vector

The pipeline runs in two passes:
- **Pass 1:** accumulate COM RDF `g(r)`.
- **Pass 2:** estimate lattice constant `a` from the **first peak** of averaged `g(r)` (unless overridden), then compute `gT(r)` using that `a`.

You can optionally provide a manual `A_LATTICE_OVERRIDE` as the final CLI argument.

## Build

```bash
make
```

## Run

```bash
./translational_order DATA_DIR START_IDX END_IDX OUT_DIR LBOND DR USE_PBC BOX_X BOX_Y [A_LATTICE_OVERRIDE]
```

Examples:

```bash
# auto-estimate a from first RDF peak
./translational_order ./data/ 1000 2000 ./out/ 1.5 0.5 1 180.0 180.0

# manual override for a
./translational_order ./data/ 1000 2000 ./out/ 1.5 0.5 1 180.0 180.0 1.12
```
