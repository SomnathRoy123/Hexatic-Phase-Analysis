#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glob.h>

#include "utils.h"
#include "io.h"
#include "clusters.h"
#include "com.h"
#include "delaunay.h"
#include "psi6.h"
#include "graccum.h"
#include "gtaccum.h"

static int extract_time_index(const char *path){
    const char *p = strrchr(path, '/');
    p = p ? p + 1 : path;
    int idx = -1;
    return (sscanf(p, "time_%d.dat", &idx) == 1) ? idx : -1;
}

static int cmp_paths_by_time(const void *a, const void *b){
    const char *pa = *(const char * const *)a;
    const char *pb = *(const char * const *)b;
    int ia = extract_time_index(pa);
    int ib = extract_time_index(pb);
    if(ia < ib) return -1;
    if(ia > ib) return 1;
    return strcmp(pa, pb);
}

static void ensure_output_dir(const char *outdir){
    char cmd[4096];
    snprintf(cmd, sizeof(cmd), "mkdir -p '%s'", outdir);
    (void)system(cmd);
}

static void usage(const char *prog){
    fprintf(stderr,
            "Usage: %s DATA_DIR START_IDX END_IDX OUT_DIR LBOND DR A_LATTICE USE_PBC BOX_X BOX_Y\n"
            "Example: %s ./snapshots/ 1000 2000 ./out/ 1.5 0.5 1.12 1 180.0 180.0\n",
            prog, prog);
}

int main(int argc, char **argv){
    if(argc < 11){
        usage(argv[0]);
        return 1;
    }

    const char *data_dir = argv[1];
    int start_idx = atoi(argv[2]);
    int end_idx = atoi(argv[3]);
    const char *out_dir = argv[4];
    double lbond = atof(argv[5]);
    double dr = atof(argv[6]);
    double a_lattice = atof(argv[7]);
    int use_pbc_flag = atoi(argv[8]);
    double box_x = atof(argv[9]);
    double box_y = atof(argv[10]);

    if(start_idx > end_idx){
        fprintf(stderr, "start index (%d) > end index (%d)\n", start_idx, end_idx);
        return 1;
    }
    if(dr <= 0.0 || lbond <= 0.0 || a_lattice <= 0.0){
        fprintf(stderr, "lbond, dr, a_lattice must all be > 0\n");
        return 1;
    }
    if(use_pbc_flag && (box_x <= 0.0 || box_y <= 0.0)){
        fprintf(stderr, "box_x and box_y must be > 0 when USE_PBC=1\n");
        return 1;
    }

    ensure_output_dir(out_dir);

    char pattern[4096];
    snprintf(pattern, sizeof(pattern), "%stime_*.dat", data_dir);

    glob_t g;
    memset(&g, 0, sizeof(g));
    if(glob(pattern, 0, NULL, &g) != 0){
        fprintf(stderr, "No files matched: %s\n", pattern);
        globfree(&g);
        return 1;
    }

    char **paths = (char**)malloc(g.gl_pathc * sizeof(char*));
    if(!paths){
        fprintf(stderr, "OOM for path list\n");
        globfree(&g);
        return 1;
    }

    size_t nsel = 0;
    for(size_t i = 0; i < g.gl_pathc; ++i){
        int ti = extract_time_index(g.gl_pathv[i]);
        if(ti >= start_idx && ti <= end_idx){
            paths[nsel] = strdup(g.gl_pathv[i]);
            if(paths[nsel]) nsel++;
        }
    }
    globfree(&g);

    if(nsel == 0){
        fprintf(stderr, "No files in range [%d, %d]\n", start_idx, end_idx);
        free(paths);
        return 1;
    }

    qsort(paths, nsel, sizeof(char*), cmp_paths_by_time);

    GrAccum *Agr = graccum_create(dr);
    GtAccum *Agt = gtaccum_create(dr);
    if(!Agr || !Agt){
        fprintf(stderr, "Accumulator creation failed\n");
        graccum_free(Agr);
        gtaccum_free(Agt);
        for(size_t i = 0; i < nsel; ++i) free(paths[i]);
        free(paths);
        return 1;
    }

    for(size_t ip = 0; ip < nsel; ++ip){
        Vec2Array pos;
        v2a_init(&pos);
        if(!read_snapshot_xy(paths[ip], &pos) || pos.n == 0){
            v2a_free(&pos);
            free(paths[ip]);
            continue;
        }

        int nclusters = 0;
        int *cluster_id = find_clusters_from_vec2array(&pos, lbond, use_pbc_flag ? 1 : 0, box_x, box_y, &nclusters);
        if(!cluster_id || nclusters < 2){
            free(cluster_id);
            v2a_free(&pos);
            free(paths[ip]);
            continue;
        }

        IntArray *clusters = make_clusters_from_ids(cluster_id, (int)pos.n, nclusters);
        if(!clusters){
            free(cluster_id);
            v2a_free(&pos);
            free(paths[ip]);
            continue;
        }

        Vec2Array coms;
        v2a_init(&coms);
        if(compute_cluster_coms(&pos, clusters, nclusters, use_pbc_flag ? 1 : 0, box_x, box_y, &coms) != 0 || coms.n < 2){
            v2a_free(&coms);
            for(int k = 0; k < nclusters; ++k) ia_free(&clusters[k]);
            free(clusters);
            free(cluster_id);
            v2a_free(&pos);
            free(paths[ip]);
            continue;
        }

        int M = 0;
        IntArray *neighbors = triangulate_get_neighbors(&coms, use_pbc_flag ? 1 : 0, box_x, box_y, &M);
        if(!neighbors || M != (int)coms.n){
            neighbors_free(neighbors, M);
            v2a_free(&coms);
            for(int k = 0; k < nclusters; ++k) ia_free(&clusters[k]);
            free(clusters);
            free(cluster_id);
            v2a_free(&pos);
            free(paths[ip]);
            continue;
        }

        Complex *psi6 = compute_psi6_from_neighbors(&coms, neighbors, use_pbc_flag ? 1 : 0, box_x, box_y);
        if(!psi6){
            neighbors_free(neighbors, M);
            v2a_free(&coms);
            for(int k = 0; k < nclusters; ++k) ia_free(&clusters[k]);
            free(clusters);
            free(cluster_id);
            v2a_free(&pos);
            free(paths[ip]);
            continue;
        }

        graccum_accumulate(Agr, &coms, use_pbc_flag ? 1 : 0, box_x, box_y);
        double theta_G = compute_global_orientation_angle(psi6, M);
        gtaccum_accumulate(Agt, &coms, theta_G, a_lattice, use_pbc_flag ? 1 : 0, box_x, box_y);

        free(psi6);
        neighbors_free(neighbors, M);
        v2a_free(&coms);
        for(int k = 0; k < nclusters; ++k) ia_free(&clusters[k]);
        free(clusters);
        free(cluster_id);
        v2a_free(&pos);
        free(paths[ip]);
    }

    free(paths);

    char out_gr[4096], out_gt[4096];
    snprintf(out_gr, sizeof(out_gr), "%s/gr_avg_time_%d_%d.dat", out_dir, start_idx, end_idx);
    snprintf(out_gt, sizeof(out_gt), "%s/gt_avg_time_%d_%d.dat", out_dir, start_idx, end_idx);

    int rc = 0;
    if(graccum_write(Agr, out_gr, start_idx, end_idx, use_pbc_flag ? 1 : 0, box_x, box_y) != 0) rc = 1;
    if(gtaccum_write(Agt, out_gt, start_idx, end_idx, a_lattice, use_pbc_flag ? 1 : 0, box_x, box_y) != 0) rc = 1;

    graccum_free(Agr);
    gtaccum_free(Agt);

    return rc;
}
