#ifndef CLUSTERS_H
#define CLUSTERS_H

#include <stdbool.h>
#include "utils.h"   /* provides Vec2, Vec2Array, IntArray, mic_delta, ia_* and v2a_* APIs */

/*
 * find_clusters_from_vec2array
 *
 * Inputs:
 *   - pos: pointer to a Vec2Array containing N = pos->n particle positions
 *   - lbond: distance cutoff (inclusive) to consider two particles connected
 *   - use_pbc: if true, minimum-image is applied using box_x, box_y
 *   - box_x, box_y: box dimensions when use_pbc is true (must be > 0)
 *
 * Outputs:
 *   - returns a malloc'd int array cluster_id of length N,
 *       cluster_id[i] in [0 .. nclusters-1]
 *   - sets *out_nclusters to the number of clusters found
 *
 * Ownership:
 *   - Caller must free() the returned array when done.
 */
int *find_clusters_from_vec2array(const Vec2Array *pos,
                                  double lbond,
                                  bool use_pbc,
                                  double box_x,
                                  double box_y,
                                  int *out_nclusters);

/*
 * make_clusters_from_ids
 *
 * Inputs:
 *   - cluster_id: array length N mapping particle -> cluster index
 *   - N: number of particles
 *   - nclusters: number of clusters (max cluster_id + 1)
 *
 * Returns:
 *   - malloc'd pointer to an array of IntArray of length nclusters.
 *     For cluster k, arr[k].data[0..arr[k].n-1] are member particle indices.
 *
 * Ownership:
 *   - Caller must call ia_free(&arr[k]) for k=0..nclusters-1, then free(arr).
 */
IntArray *make_clusters_from_ids(const int *cluster_id, int N, int nclusters);

#endif /* CLUSTERS_H */
