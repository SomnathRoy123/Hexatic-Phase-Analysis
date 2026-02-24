#ifndef DELAUNAY_H
#define DELAUNAY_H

#include "utils.h"
#include <stdbool.h>

/*
 * triangulate_get_neighbors
 *
 * Given `points` (Vec2Array, M points), compute Delaunay neighbors for the M original points.
 *
 * If use_pbc is true, 8 periodic image copies are created so edges crossing the box
 * boundaries are found. The function returns an array of IntArray of length M:
 *   neighbors[i] contains the unique neighbor indices (in [0..M-1]) of point i.
 *
 * On success: returns pointer to malloc'd IntArray array (length M). Caller must call
 *    neighbors_free(neighbors, M);
 * On failure: returns NULL.
 *
 * Note: this function uses Triangle (triangulate). It will free Triangle-allocated
 * memory via trifree() and free() for memory it allocates.
 */
IntArray *triangulate_get_neighbors(const Vec2Array *points,
                                    bool use_pbc,
                                    double box_x, double box_y,
                                    int *out_M);

/* Free neighbors array returned by triangulate_get_neighbors */
void neighbors_free(IntArray *neighbors, int M);

#endif /* DELAUNAY_H */
