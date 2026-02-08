
/*
 * clusters.c
 *
 * Implements cluster finding (union-find) on 2D point sets with optional PBC.
 *
 * Compile: include utils.c/utils.h in your build and compile with -std=c99 -O2 -Wall
 */

#include "clusters.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* -------------------- Internal Union-Find -------------------- */
typedef struct {
    int *parent;
    int *rank;
    int  n;
} UF;

static void uf_init(UF *uf, int n){
    uf->n = n;
    uf->parent = (int*)malloc(n * sizeof(int));
    uf->rank   = (int*)calloc(n, sizeof(int));
    if(!uf->parent || !uf->rank){
        fprintf(stderr, "UF: out of memory\n");
        exit(1);
    }
    for(int i=0;i<n;i++) uf->parent[i] = i;
}

static void uf_free(UF *uf){
    if(uf->parent){
        free(uf->parent);
    } 
    if(uf->rank){
        free(uf->rank);
    }
    uf->parent = NULL;
    uf->rank = NULL;
    uf->n = 0;
}

static int uf_root(UF *uf, int i){
    /* path compression */
    int r = i;
    while(uf->parent[r] != r){
         r = uf->parent[r];
    }

    while(i != r){
        int p = uf->parent[i];
        uf->parent[i] = r;
        i = p;
    }
    return r;
}

static void uf_union(UF *uf, int a, int b){
    int ra = uf_root(uf, a);
    int rb = uf_root(uf, b);
    if(ra == rb) return;
    if(uf->rank[ra] < uf->rank[rb]){
        uf->parent[ra] = rb;
    }
    else if(uf->rank[ra] > uf->rank[rb]){
        uf->parent[rb] = ra;
    }
    else { 
        uf->parent[rb] = ra; 
        uf->rank[ra]++; 
    }
}

/* ------------------- Public API ------------------- */

int *find_clusters_from_vec2array(const Vec2Array *pos,
                                  double lbond,
                                  bool use_pbc,
                                  double box_x,
                                  double box_y,
                                  int *out_nclusters)
{
    if(!pos || out_nclusters == NULL){
        fprintf(stderr, "find_clusters: invalid arguments\n");
        return NULL;
    }

    int N = (int)pos->n;
    *out_nclusters = 0;

    if(N <= 0){
        return NULL;
    }

    if(use_pbc && (box_x <= 0.0 || box_y <= 0.0)){
        fprintf(stderr, "find_clusters: use_pbc true but box_x/box_y not positive\n");
        return NULL;
    }

    UF uf;
    uf_init(&uf, N);
    const double lb2 = lbond * lbond;

    /* naive O(N^2) pair scan. Replace by cell-list if N large. */
    for(int i=0;i<N-1;i++){
        for(int j=i+1;j<N;j++){
            double dx = pos->data[j].x - pos->data[i].x;
            double dy = pos->data[j].y - pos->data[i].y;
            if(use_pbc){
                dx = mic_delta(dx, box_x);
                dy = mic_delta(dy, box_y);
            }
            double d2 = dx*dx + dy*dy;
            if(d2 <= lb2){
                uf_union(&uf, i, j);
            }
        }
    }

    /* compute root for each particle */
    int *root = (int*)malloc(N * sizeof(int));
    if(!root){
        fprintf(stderr,"find_clusters: OOM\n"); 
        uf_free(&uf); 
        return NULL; 
    }
    for(int i=0;i<N;i++){ 
        root[i] = uf_root(&uf, i); 
    }

    /* map distinct roots to compact cluster indices */
    int *map = (int*)malloc(N * sizeof(int));
    if(!map){
        fprintf(stderr,"find_clusters: OOM\n"); 
        free(root); 
        uf_free(&uf); 
        return NULL; 
    }

    for(int i=0;i<N;i++){
        map[i] = -1;
    }

    int nclusters = 0;
    for(int i=0;i<N;i++){
        int r = root[i];
        if(map[r] == -1) map[r] = nclusters++;
    }

    /* create cluster_id array (caller frees) */
    int *cluster_id = (int*)malloc(N * sizeof(int));
    if(!cluster_id){ fprintf(stderr,"find_clusters: OOM\n"); free(root); free(map); uf_free(&uf); return NULL; }
    for(int i=0;i<N;i++){
        cluster_id[i] = map[root[i]];
    }

    free(root);
    free(map);
    uf_free(&uf);

    *out_nclusters = nclusters;
    return cluster_id;
}

IntArray *make_clusters_from_ids(const int *cluster_id, int N, int nclusters){
    if(nclusters <= 0) return NULL;
    if(!cluster_id || N <= 0) return NULL;

    IntArray *clusters = (IntArray*)malloc(nclusters * sizeof(IntArray));
    if(!clusters){ fprintf(stderr,"make_clusters_from_ids: OOM\n"); return NULL; }
    for(int k=0;k<nclusters;k++) ia_init(&clusters[k]);

    for(int i=0;i<N;i++){
        int cid = cluster_id[i];
        if(cid < 0 || cid >= nclusters){
            fprintf(stderr,"make_clusters_from_ids: warning: invalid cid %d at i=%d (skipped)\n", cid, i);
            continue;
        }
        ia_push(&clusters[cid], i);
    }
    return clusters;
}
