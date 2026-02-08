/*
 * main.c
 *
 * Orchestrates the pipeline:
 *   for each snapshot time_<idx>.dat in DATA_DIR with idx in [START_INDEX..END_INDEX]:
 *     - read positions (x y) via read_snapshot_xy()
 *     - cluster via union-find
 *     - compute cluster COMs
 *     - compute Delaunay neighbors (with PBC image trick)
 *     - compute psi6 at COMs
 *     - accumulate g6(r)
 *   finally write averaged g6 to OUTPUT_DIR/g6_avg_time_<start>_<end>.dat
 *
 * Assumes the following modules exist and are linked:
 *   - utils.{c,h}
 *   - io.{c,h}         (read_snapshot_xy)
 *   - clusters.{c,h}
 *   - com.{c,h}
 *   - delaunay.{c,h}
 *   - psi6.{c,h}
 *   - g6accum.{c,h}
 *
 * Compile: see Makefile in project root (link everything together).
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glob.h>
#include <errno.h>

#include "utils.h"
#include "clusters.h"
#include "com.h"
#include "delaunay.h"
#include "psi6.h"
#include "g6accum.h"
#include "io.h"   /* must provide: int read_snapshot_xy(const char *path, Vec2Array *pos); */

/* ----------------------- DEFAULT CONFIG (can be moved to params.h) ----------------------- */
static const char *DEFAULT_DATA_DIR   = "/home/somnath2/Codes/Trial_102/TA160_R3.9/evo/";   /* trailing slash allowed */
static const char *DEFAULT_OUTPUT_DIR = "/home/somnath2/Codes/cluster_analysis/Hexatic Order spatial Code/data_out/TA160_R3.9/";
static const int   DEFAULT_START_IDX  = 1000;
static const int   DEFAULT_END_IDX    = 2000;
static const double DEFAULT_LBOND     = 1.5;
static const double DEFAULT_DR        = 0.5;
static const int    VERBOSITY         = 1;
/* ----------------------------------------------------------------------------------------- */

/* Helper: extract time index from filename like .../time_<num>.dat */
static int extract_time_index(const char *path){
    const char *p = strrchr(path, '/');
    if(p) p++; else p = path;
    int idx = -1;
    if(sscanf(p, "time_%d.dat", &idx) == 1) return idx;
    return -1;
}

static int cmp_paths_by_time(const void *a, const void *b){
    const char *pa = *(const char * const *)a;
    const char *pb = *(const char * const *)b;
    int ia = extract_time_index(pa);
    int ib = extract_time_index(pb);
    if(ia < ib) return -1;
    if(ia > ib) return 1;
    /* fallback lexical */
    return strcmp(pa, pb);
}

/* Ensure output dir exists (simple - uses system mkdir -p). Could be replaced by portable code. */
static void ensure_output_dir(const char *outdir){
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "mkdir -p '%s'", outdir);
    int rc = system(cmd);
    (void)rc;
}

/* Print usage */
static void usage(const char *prog){
    fprintf(stderr,
        "Usage: %s DATA_DIR START_INDEX END_INDEX OUTPUT_DIR [LBOND] [DR] [USE_PBC] [BOX_X] [BOX_Y]\n\n"
        "Example:\n"
        "  %s ./data/ 1000 1200 ./out/ 1.5 0.5 1 180.0 180.0\n\n"
        "If optional args omitted, defaults are used.\n",
        prog, prog);
}

/* ------------------------------- main ---------------------------------- */
int main(int argc, char **argv){
    const char *data_dir = DEFAULT_DATA_DIR;
    const char *out_dir  = DEFAULT_OUTPUT_DIR;
    int start_idx = DEFAULT_START_IDX;
    int end_idx   = DEFAULT_END_IDX;
    double lbond  = DEFAULT_LBOND;
    double dr     = DEFAULT_DR;
    int use_pbc_flag = 1;
    double box_x = 180.0, box_y = 180.0;

    if(argc < 5){
        if(argc == 1){
            fprintf(stderr, "No args supplied — using defaults. To see usage, run with -h\n");
        } else {
            usage(argv[0]);
            return 1;
        }
    } else {
        data_dir = argv[1];
        start_idx = atoi(argv[2]);
        end_idx   = atoi(argv[3]);
        out_dir   = argv[4];
        if(argc >= 6) lbond = atof(argv[5]);
        if(argc >= 7) dr    = atof(argv[6]);
        if(argc >= 8) use_pbc_flag = atoi(argv[7]);
        if(argc >= 10) {
            box_x = atof(argv[8]);
            box_y = atof(argv[9]);
        }
    }

    if(start_idx > end_idx){
        fprintf(stderr, "start index (%d) > end index (%d)\n", start_idx, end_idx);
        return 1;
    }

    ensure_output_dir(out_dir);

    /* Build glob pattern */
    char pattern[4096];
    snprintf(pattern, sizeof(pattern), "%stime_*.dat", data_dir);

    glob_t g;
    memset(&g, 0, sizeof(g));
    int glob_rc = glob(pattern, 0, NULL, &g);
    if(glob_rc != 0){
        fprintf(stderr, "No files match pattern: %s  (glob rc=%d)\n", pattern, glob_rc);
        globfree(&g);
        return 1;
    }

    /* Filter paths by index range */
    char **paths = (char**)malloc(g.gl_pathc * sizeof(char*));
    if(!paths){ fprintf(stderr,"OOM\n"); globfree(&g); return 1; }
    size_t nsel = 0;
    for(size_t i=0;i<g.gl_pathc;i++){
        int ti = extract_time_index(g.gl_pathv[i]);
        if(ti < 0) continue;
        if(ti >= start_idx && ti <= end_idx){
            paths[nsel++] = strdup(g.gl_pathv[i]);
        }
    }
    globfree(&g);

    if(nsel == 0){
        fprintf(stderr, "No files found in range %d..%d\n", start_idx, end_idx);
        free(paths);
        return 1;
    }

    qsort(paths, nsel, sizeof(char*), cmp_paths_by_time);
    if(VERBOSITY) printf("Found %zu files in range [%d, %d]\n", nsel, start_idx, end_idx);

    /* Create accumulator */
    G6Accum *A = g6accum_create(dr);
    if(!A){ fprintf(stderr,"Failed to create g6 accumulator\n"); for(size_t i=0;i<nsel;i++) free(paths[i]); free(paths); return 1; }

    /* Process snapshots */
    for(size_t ip=0; ip<nsel; ip++){
        const char *path = paths[ip];
        int tindex = extract_time_index(path);
        if(VERBOSITY) printf("[%zu/%zu] Processing %s (t=%d)\n", ip+1, nsel, path, tindex);






                /* 1) Read snapshot positions (expects io.c to implement read_snapshot_xy) */
        Vec2Array pos;
        v2a_init(&pos);
        if(!read_snapshot_xy(path, &pos)){
            fprintf(stderr, "  ! failed to read %s (skipping)\n", path);
            v2a_free(&pos);
            free((void*)path);
            continue;
        }
        printf("  read %zu particles\n", pos.n);   // <-- ADD THIS

        if(pos.n == 0){
            fprintf(stderr, "  ! empty snapshot %s (skipping)\n", path);
            v2a_free(&pos);
            free((void*)path);
            continue;
        }

        /* 2) Clustering (union-find) */
        printf("  entering clustering\n");         // <-- ADD THIS
        int nclusters = 0;
        int *cluster_id = find_clusters_from_vec2array(&pos, lbond, use_pbc_flag ? 1 : 0, box_x, box_y, &nclusters);
        printf("  clustering done, nclusters = %d\n", nclusters);   // <-- ADD THIS


        // /* 1) Read snapshot positions (expects io.c to implement read_snapshot_xy) */
        // Vec2Array pos;
        // v2a_init(&pos);
        // if(!read_snapshot_xy(path, &pos)){
        //     fprintf(stderr, "  ! failed to read %s (skipping)\n", path);
        //     v2a_free(&pos);
        //     free((void*)path);
        //     continue;
        // }
        // if(pos.n == 0){
        //     fprintf(stderr, "  ! empty snapshot %s (skipping)\n", path);
        //     v2a_free(&pos);
        //     free((void*)path);
        //     continue;
        // }

        // /* 2) Clustering (union-find) */
        // int nclusters = 0;
        // int *cluster_id = find_clusters_from_vec2array(&pos, lbond, use_pbc_flag ? 1 : 0, box_x, box_y, &nclusters);


        /* 2) Clustering (union-find) */


        if (!cluster_id) {
            fprintf(stderr, "  ! clustering failed (null cluster_id)\n");
            v2a_free(&pos);
            free((void*)path);
            continue;
        }

        /* DEBUG: check label range */
        int max_id = -1, min_id = 1e9;
        for (int i = 0; i < (int)pos.n; i++) {
            if (cluster_id[i] < min_id) min_id = cluster_id[i];
            if (cluster_id[i] > max_id) max_id = cluster_id[i];
        }
        printf("  cluster_id range: [%d, %d]\n", min_id, max_id);
        if (min_id < 0 || max_id >= nclusters) {
            fprintf(stderr,
                    "  !! ERROR: cluster_id out of range: min=%d max=%d nclusters=%d\n",
                    min_id, max_id, nclusters);
            /* bail out so we see the message instead of segfault */
            free(cluster_id);
            v2a_free(&pos);
            free((void*)path);
            continue;
        }

        // if(!cluster_id){
        //     fprintf(stderr, "  ! clustering failed (skipping)\n");
        //     v2a_free(&pos);
        //     free((void*)path);
        //     continue;
        // }
        // if(nclusters < 2){
        //     if(VERBOSITY) printf("  ! less than 2 clusters (%d) — skipping\n", nclusters);
        //     free(cluster_id);
        //     v2a_free(&pos);
        //     free((void*)path);
        //     continue;
        // }

        /* 3) Build IntArray clusters and compute COMs */
        printf("  building clusters (make_clusters_from_ids)\n");
        IntArray *clusters = make_clusters_from_ids(cluster_id, (int)pos.n, nclusters);
        if (!clusters) {
            fprintf(stderr, "  ! make_clusters_from_ids returned NULL\n");
            free(cluster_id);
            v2a_free(&pos);
            free((void*)path);
            continue;
        }
        printf("  clusters built\n");

        
        Vec2Array coms;
        v2a_init(&coms);
        printf("  computing COMs\n");
        if(compute_cluster_coms(&pos, clusters, nclusters, use_pbc_flag ? 1 : 0, box_x, box_y, &coms) != 0){
            fprintf(stderr, "  ! compute_cluster_coms failed (skipping)\n");
            for(int k=0;k<nclusters;k++) ia_free(&clusters[k]);
            free(clusters);
            free(cluster_id);
            v2a_free(&pos);
            free((void*)path);
            continue;
        }
        printf("  COMs computed: %zu clusters\n", coms.n);
        /* 4) Delaunay neighbors (with PBC images) */
        int M;
        IntArray *neighbors = triangulate_get_neighbors(&coms, use_pbc_flag ? 1 : 0, box_x, box_y, &M);
        printf("  triangulation returned neighbors, M = %d\n", M);
        if(!neighbors || M != (int)coms.n){
            fprintf(stderr, "  ! triangulate_get_neighbors failed (skipping)\n");
            if(neighbors) neighbors_free(neighbors, M);
            v2a_free(&coms);
            for(int k=0;k<nclusters;k++) ia_free(&clusters[k]);
            free(clusters);
            free(cluster_id);
            v2a_free(&pos);
            free((void*)path);
            continue;
        }

        /* 5) psi6 */
        printf("  computing psi6\n");
        Complex *psi6 = compute_psi6_from_neighbors(&coms, neighbors, use_pbc_flag ? 1 : 0, box_x, box_y);
        if(!psi6){
            fprintf(stderr, "  ! compute_psi6 failed (skipping)\n");
            neighbors_free(neighbors, M);
            v2a_free(&coms);
            for(int k=0;k<nclusters;k++) ia_free(&clusters[k]);
            free(clusters);
            free(cluster_id);
            v2a_free(&pos);
            free((void*)path);
            continue;
        }
        printf("  psi6 computed\n");
        /* 6) accumulate g6 */
        g6accum_accumulate(A, &coms, psi6, use_pbc_flag ? 1 : 0, box_x, box_y);

        /* cleanup per-snapshot */
        free(psi6);
        neighbors_free(neighbors, M);

        v2a_free(&coms);
        for(int k=0;k<nclusters;k++) ia_free(&clusters[k]);
        free(clusters);

        free(cluster_id);
        v2a_free(&pos);

        free((void*)path);
    }

    free(paths);

    /* Write averaged file */
    char outpath[4096];
    snprintf(outpath, sizeof(outpath), "%s/g6_avg_time_%d_%d.dat", out_dir, start_idx, end_idx);
    if(g6accum_write(A, outpath, start_idx, end_idx, lbond, use_pbc_flag ? 1 : 0, box_x, box_y) != 0){
        fprintf(stderr, "Failed to write g6 average file\n");
        g6accum_free(A);
        return 1;
    }
    g6accum_free(A);

    if(VERBOSITY) printf("✓ Done. Wrote %s\n", outpath);
    return 0;
}
