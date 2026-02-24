/*
 * com.c
 *
 * Compute cluster centers-of-mass (COM) robustly with optional periodic boundary handling.
 *
 * Requires: utils.h (Vec2Array, IntArray, v2a_ia_* and mic_delta/wrap_pos)
 *
 * Ownership:
 *   - compute_cluster_coms initializes `coms` with v2a_init() and pushes entries.
 *   - Caller must call v2a_free(&coms).
 */

#include "com.h"
#include <stdlib.h>
#include <stdio.h>

/* Helper: safe wrap (delegates to wrap_pos which is in utils.c) */
static inline double safe_wrap(double x, double L){
    return wrap_pos(x, L);
}

/* Main implementation */
int compute_cluster_coms(const Vec2Array *pos,
                         const IntArray *clusters,
                         int nclusters,
                         bool use_pbc,
                         double box_x,
                         double box_y,
                         Vec2Array *coms)
{
    if(!pos || !clusters || !coms || nclusters < 0){
        fprintf(stderr, "compute_cluster_coms: invalid arguments\n");
        return 1;
    }
    if(use_pbc && (box_x <= 0.0 || box_y <= 0.0)){
        fprintf(stderr, "compute_cluster_coms: invalid box dims for PBC\n");
        return 2;
    }

    v2a_init(coms);

    const int N = (int)pos->n;
    for(int c = 0; c < nclusters; ++c){
        int nm = (int)clusters[c].n;
        if(nm <= 0){
            /* empty cluster -> push zero vector (caller decides how to handle) */
            v2a_push(coms, (Vec2){0.0, 0.0});
            continue;
        }

        /* Reference particle is the first member */
        int idx0 = clusters[c].data[0];
        if(idx0 < 0 || idx0 >= N){
            fprintf(stderr, "compute_cluster_coms: invalid member index %d in cluster %d\n", idx0, c);
            v2a_push(coms, (Vec2){0.0, 0.0});
            continue;
        }
        double x0 = pos->data[idx0].x;
        double y0 = pos->data[idx0].y;

        double sx = 0.0;
        double sy = 0.0;

        for(int k = 0; k < nm; ++k){
            int idx = clusters[c].data[k];
            if(idx < 0 || idx >= N){
                fprintf(stderr, "compute_cluster_coms: invalid member index %d in cluster %d (skipped)\n", idx, c);
                continue;
            }
            double dx = pos->data[idx].x - x0;
            double dy = pos->data[idx].y - y0;
            if(use_pbc){
                dx = mic_delta(dx, box_x);
                dy = mic_delta(dy, box_y);
            }
            sx += dx;
            sy += dy;
        }

        double cx = x0 + sx / (double)nm;
        double cy = y0 + sy / (double)nm;

        if(use_pbc){
            cx = safe_wrap(cx, box_x);
            cy = safe_wrap(cy, box_y);
        }

        v2a_push(coms, (Vec2){cx, cy});
    }

    return 0;
}

/* Convenience wrapper that builds IntArray clusters from cluster_id, uses compute_cluster_coms */
int compute_cluster_coms_from_ids(const Vec2Array *pos,
                                  const int *cluster_id, int N, int nclusters,
                                  bool use_pbc, double box_x, double box_y,
                                  Vec2Array *coms)
{
    if(!pos || !cluster_id || !coms || N < 0 || nclusters < 0){
        fprintf(stderr, "compute_cluster_coms_from_ids: invalid args\n");
        return 1;
    }

    /* build IntArray list */
    IntArray *clusters = (IntArray*)malloc((size_t)nclusters * sizeof(IntArray));
    if(!clusters){ fprintf(stderr, "compute_cluster_coms_from_ids: OOM\n"); return 2; }
    for(int k=0;k<nclusters;k++) ia_init(&clusters[k]);

    for(int i=0;i<N;i++){
        int cid = cluster_id[i];
        if(cid < 0 || cid >= nclusters){
            fprintf(stderr, "compute_cluster_coms_from_ids: invalid cid %d at i=%d (skipped)\n", cid, i);
            continue;
        }
        ia_push(&clusters[cid], i);
    }

    int rc = compute_cluster_coms(pos, clusters, nclusters, use_pbc, box_x, box_y, coms);

    for(int k=0;k<nclusters;k++) ia_free(&clusters[k]);
    free(clusters);

    return rc;
}
