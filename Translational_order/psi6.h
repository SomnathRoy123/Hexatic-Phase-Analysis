#ifndef PSI6_H
#define PSI6_H

#include "utils.h"
#include "delaunay.h"
#include <stdbool.h>

/* Complex type (small) */
typedef struct { double re, im; } Complex;

/*
 * compute_psi6_from_neighbors
 *
 * Inputs:
 *   - coms: Vec2Array containing positions of original points (M points)
 *   - neighbors: IntArray array (length M) where neighbors[i] lists neighbor indices (0..M-1)
 *   - use_pbc, box_x, box_y: for angle calculation (apply minimum image if use_pbc true)
 *
 * Returns:
 *   - malloc'd Complex array of length M where psi6[i] holds the normalized local psi6.
 *     Caller must free() the result.
 *   - On failure: returns NULL.
 *
 * Behavior:
 *   - For a point with zero neighbors, psi6 is {0,0}.
 *   - Neighbor lists are assumed already deduplicated.
 */
Complex *compute_psi6_from_neighbors(const Vec2Array *coms,
                                     const IntArray *neighbors,
                                     bool use_pbc,
                                     double box_x, double box_y);

double compute_global_orientation_angle(const Complex *psi6, int M);

#endif /* PSI6_H */
