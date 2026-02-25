/*
 * delaunay.c
 *
 * Triangle wrapper that returns deduplicated neighbor lists for original points.
 *
 * Requires triangle.h in include path and Triangle library (or triangle.c) linked.
 *
 * Important: free() memory you allocate, and use trifree() only on Triangle-allocated outputs.
 */

#include "delaunay.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define REAL double
#include "triangle.h"  /* Triangle API; ensure this is available at compile time */

/* Utility: check for neighbor presence (linear) */
static inline int neighbor_has(const IntArray *nbrs, int v){
    for(size_t k=0;k<nbrs->n;k++){
        if(nbrs->data[k] == v) return 1;
    }
    return 0;
}

IntArray *triangulate_get_neighbors(const Vec2Array *points,
                                    bool use_pbc,
                                    double box_x, double box_y,
                                    int *out_M)
{
    if(!points || !out_M){ 
        fprintf(stderr,"triangulate: invalid args\n");
        return NULL;
     }
    const int M = (int)points->n;
    *out_M = 0;
    if(M <= 0){
        return NULL;
     }

    /* build Triangle input structure (points + images if PBC) */
    struct triangulateio in, out;
    memset(&in, 0, sizeof(in));
    memset(&out, 0, sizeof(out));

    int total_points = M;
    if(use_pbc) {
        if(box_x <= 0.0 || box_y <= 0.0){ fprintf(stderr,"triangulate: invalid box dims\n"); return NULL; }
        total_points = M * 9; /* original + 8 images */
    }

    in.numberofpoints = total_points;
    in.pointlist = (REAL*)malloc((size_t)total_points * 2 * sizeof(REAL));
    if(!in.pointlist){ fprintf(stderr,"triangulate: OOM\n"); return NULL; }

    /* fill originals */
    for(int i=0;i<M;i++){
        in.pointlist[i*2 + 0] = points->data[i].x;
        in.pointlist[i*2 + 1] = points->data[i].y;
    }

    /* fill images */
    if(use_pbc){
        int idx = M;
        int shifts[8][2] = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,1},{1,-1},{1,0},{1,1}};
        for(int s=0;s<8;s++){
            double sx = shifts[s][0] * box_x;
            double sy = shifts[s][1] * box_y;
            for(int i=0;i<M;i++){
                in.pointlist[idx*2 + 0] = points->data[i].x + sx;
                in.pointlist[idx*2 + 1] = points->data[i].y + sy;
                idx++;
            }
        }
    }

    in.numberofpointattributes = 0;
    in.pointmarkerlist = NULL;
    in.numberofsegments = 0;
    in.segmentlist = NULL;

    /* Run Triangle: z = zero-based indexing, P = no new points, E = output edges, Q = quiet */
    char switches[] = "zeQ";
    triangulate(switches, &in, &out, NULL);

    if(out.numberofedges <= 0){
        fprintf(stderr,"triangulate: no edges produced\n");
        free(in.pointlist);
        if(out.pointlist) trifree(out.pointlist);
        if(out.edgelist) trifree(out.edgelist);
        if(out.trianglelist) trifree(out.trianglelist);
        return NULL;
    }

    /* Allocate neighbor arrays for original M points */
    IntArray *neighbors = (IntArray*)malloc(sizeof(IntArray) * (size_t)M);
    if(!neighbors){ fprintf(stderr,"triangulate: OOM neighbors\n"); free(in.pointlist); if(out.pointlist) trifree(out.pointlist); if(out.edgelist) trifree(out.edgelist); if(out.trianglelist) trifree(out.trianglelist); return NULL; }
    for(int i=0;i<M;i++) ia_init(&neighbors[i]);

    /*
     * Scan edges and add unique neighbor pairs.
     *
     * For PBC runs we triangulate a 3x3 tiled point set (original + 8 images).
     * We must only keep edges touching the central (original) tile; otherwise
     * edges entirely between image tiles can fold back (via %M) into artificial
     * long-range neighbors in the original box.
     */
    for(int e=0; e < out.numberofedges; e++){
        int p1 = out.edgelist[e*2 + 0];
        int p2 = out.edgelist[e*2 + 1];

        if(use_pbc){
            int p1_is_central = (p1 < M);
            int p2_is_central = (p2 < M);
            if(!p1_is_central && !p2_is_central) continue;
        }

        /* Map to original index [0..M-1] */
        int o1 = p1 % M;
        int o2 = p2 % M;

        /* Ignore self-edge */
        if(o1 == o2) continue;

        if(use_pbc){
            /*
             * Keep only edges consistent with the minimum-image displacement
             * between the mapped originals. This rejects central-image links
             * to a non-nearest periodic copy that can appear in the tiled mesh.
             */
            double ex = in.pointlist[p2*2 + 0] - in.pointlist[p1*2 + 0];
            double ey = in.pointlist[p2*2 + 1] - in.pointlist[p1*2 + 1];

            double mx = points->data[o2].x - points->data[o1].x;
            double my = points->data[o2].y - points->data[o1].y;
            mx = mic_delta(mx, box_x);
            my = mic_delta(my, box_y);

            const double eps = 1e-9;
            if(fabs(ex - mx) > eps || fabs(ey - my) > eps) continue;
        }

        /* Add o2 to neighbors of o1 if not present */
        if(!neighbor_has(&neighbors[o1], o2)) ia_push(&neighbors[o1], o2);

        /* Add o1 to neighbors of o2 if not present */
        if(!neighbor_has(&neighbors[o2], o1)) ia_push(&neighbors[o2], o1);
    }

    /* Cleanup Triangle memory. in.pointlist was malloc'ed by us -> free() it.
       out.* must be freed with trifree() if non-NULL. */
    free(in.pointlist);
    if(out.pointlist) trifree(out.pointlist);
    if(out.edgelist) trifree(out.edgelist);
    if(out.trianglelist) trifree(out.trianglelist);

    *out_M = M;
    return neighbors;
}

void neighbors_free(IntArray *neighbors, int M){
    if(!neighbors) return;
    for(int i=0;i<M;i++) ia_free(&neighbors[i]);
    free(neighbors);
}
