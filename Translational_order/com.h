#ifndef COM_H
#define COM_H

#include <stdbool.h>
#include "utils.h"   /* Vec2Array, IntArray, v2a_* and ia_* APIs */

/*
 * compute_cluster_coms
 *
 * Inputs:
 *   - pos: Vec2Array of particle positions (length N = pos->n)
 *   - clusters: array of IntArray (length nclusters), each lists member indices
 *   - nclusters: number of clusters
 *   - use_pbc: whether to apply minimum-image when computing relative displacements
 *   - box_x, box_y: required when use_pbc is true (must be >0)
 *
 * Output:
 *   - coms: an output Vec2Array (UNINITIALIZED by caller). This function will call
 *           v2a_init(&coms) and push `nclusters` entries. Caller must v2a_free(&coms).
 *
 * Returns:
 *   - 0 on success, non-zero on error (e.g., invalid args).
 */
int compute_cluster_coms(const Vec2Array *pos,
                         const IntArray *clusters,
                         int nclusters,
                         bool use_pbc,
                         double box_x,
                         double box_y,
                         Vec2Array *coms);

/*
 * compute_cluster_coms_from_ids
 *
 * Convenience wrapper:
 * Inputs:
 *   - pos: Vec2Array of particle positions
 *   - cluster_id: int array length N mapping particle -> cluster index (0..nclusters-1)
 *   - N: number of particles
 *   - nclusters: number of clusters
 *   - use_pbc, box_x, box_y: PBC settings
 *
 * Output:
 *   - coms: output Vec2Array (v2a_init called inside); caller must v2a_free(&coms)
 *
 * Returns:
 *   - 0 on success, non-zero on error.
 */
int compute_cluster_coms_from_ids(const Vec2Array *pos,
                                  const int *cluster_id, int N, int nclusters,
                                  bool use_pbc, double box_x, double box_y,
                                  Vec2Array *coms);

#endif /* COM_H */
