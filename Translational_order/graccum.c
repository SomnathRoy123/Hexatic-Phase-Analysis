#include "graccum.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    double r_center;
    double shell_area_sum;
    double ideal_pairs_sum;
    long pair_count;
} GrBin;

struct GrAccum {
    GrBin *bins;
    int nbins;
    double dr;
    long nframes;
};

static double shell_area(double rin, double rout){
    return M_PI * (rout * rout - rin * rin);
}

GrAccum *graccum_create(double dr){
    if(dr <= 0.0){
        fprintf(stderr, "graccum_create: dr must be > 0\n");
        return NULL;
    }
    GrAccum *A = (GrAccum*)malloc(sizeof(GrAccum));
    if(!A){
        fprintf(stderr, "graccum_create: OOM\n");
        return NULL;
    }
    A->bins = NULL;
    A->nbins = 0;
    A->dr = dr;
    A->nframes = 0;
    return A;
}

void graccum_free(GrAccum *A){
    if(!A) return;
    free(A->bins);
    free(A);
}

static void ensure_bins(GrAccum *A, int bmax){
    if(bmax < 0 || bmax < A->nbins) return;
    int new_n = bmax + 1;
    GrBin *nb = (GrBin*)realloc(A->bins, (size_t)new_n * sizeof(GrBin));
    if(!nb){
        fprintf(stderr, "graccum_ensure_bins: OOM\n");
        exit(1);
    }
    for(int b = A->nbins; b < new_n; ++b){
        nb[b].r_center = (b + 0.5) * A->dr;
        nb[b].shell_area_sum = 0.0;
        nb[b].ideal_pairs_sum = 0.0;
        nb[b].pair_count = 0;
    }
    A->bins = nb;
    A->nbins = new_n;
}

void graccum_accumulate(GrAccum *A,
                        const Vec2Array *coms,
                        bool use_pbc,
                        double box_x,
                        double box_y){
    if(!A || !coms) return;
    int M = (int)coms->n;
    if(M < 2) return;

    if(use_pbc && (box_x <= 0.0 || box_y <= 0.0)){
        fprintf(stderr, "graccum_accumulate: invalid PBC box\n");
        return;
    }

    double rmax2 = 0.0;
    for(int i = 0; i < M - 1; ++i){
        for(int j = i + 1; j < M; ++j){
            double dx = coms->data[j].x - coms->data[i].x;
            double dy = coms->data[j].y - coms->data[i].y;
            if(use_pbc){
                dx = mic_delta(dx, box_x);
                dy = mic_delta(dy, box_y);
            }
            double r2 = dx*dx + dy*dy;
            if(r2 > rmax2) rmax2 = r2;
        }
    }

    int bmax = (int)floor(sqrt(rmax2) / A->dr);
    ensure_bins(A, bmax);

    for(int i = 0; i < M - 1; ++i){
        for(int j = i + 1; j < M; ++j){
            double dx = coms->data[j].x - coms->data[i].x;
            double dy = coms->data[j].y - coms->data[i].y;
            if(use_pbc){
                dx = mic_delta(dx, box_x);
                dy = mic_delta(dy, box_y);
            }
            double r = sqrt(dx*dx + dy*dy);
            int b = (int)floor(r / A->dr);
            if(b >= 0 && b < A->nbins) A->bins[b].pair_count += 1;
        }
    }

    double area = use_pbc ? box_x * box_y : 0.0;
    double rho = (area > 0.0) ? ((double)M / area) : 0.0;

    for(int b = 0; b < A->nbins; ++b){
        double rin = b * A->dr;
        double rout = (b + 1) * A->dr;
        double da = shell_area(rin, rout);
        A->bins[b].shell_area_sum += da;
        if(rho > 0.0){
            A->bins[b].ideal_pairs_sum += 0.5 * (double)M * rho * da;
        }
    }

    A->nframes += 1;
}


static double gr_value(const GrBin *bin){
    if(!bin) return 0.0;
    if(bin->ideal_pairs_sum <= 0.0) return 0.0;
    return (double)bin->pair_count / bin->ideal_pairs_sum;
}

double graccum_estimate_first_peak_a(const GrAccum *A){
    if(!A || !A->bins || A->nbins < 3) return -1.0;

    int first_local_peak = -1;
    int global_peak = -1;
    double global_peak_val = -1.0;

    for(int b = 1; b < A->nbins - 1; ++b){
        double gm = gr_value(&A->bins[b - 1]);
        double g0 = gr_value(&A->bins[b]);
        double gp = gr_value(&A->bins[b + 1]);

        if(g0 > global_peak_val){
            global_peak_val = g0;
            global_peak = b;
        }

        if(first_local_peak < 0 && g0 > gm && g0 >= gp && g0 > 1.0){
            first_local_peak = b;
        }
    }

    int pick = (first_local_peak >= 0) ? first_local_peak : global_peak;
    if(pick < 0 || pick >= A->nbins) return -1.0;

    return A->bins[pick].r_center;
}

int graccum_write(GrAccum *A,
                  const char *outpath,
                  int t0,
                  int t1,
                  bool use_pbc,
                  double box_x,
                  double box_y){
    if(!A || !outpath) return 1;

    FILE *f = fopen(outpath, "w");
    if(!f){
        fprintf(stderr, "graccum_write: cannot open %s: %s\n", outpath, strerror(errno));
        return 2;
    }

    fprintf(f, "# g(r) average over time_%d..time_%d\n", t0, t1);
    fprintf(f, "# columns: r_center pair_count shell_area pair_density ideal_pairs g_r\n");
    fprintf(f, "# dr=%.8g use_pbc=%s frames=%ld\n", A->dr, use_pbc ? "true" : "false", A->nframes);
    if(use_pbc) fprintf(f, "# box: %.8g %.8g\n", box_x, box_y);

    for(int b = 0; b < A->nbins; ++b){
        double pairs = (double)A->bins[b].pair_count;
        double sa = A->bins[b].shell_area_sum;
        double pd = (sa > 0.0) ? pairs / sa : 0.0;
        double ideal = A->bins[b].ideal_pairs_sum;
        double gr = (ideal > 0.0) ? pairs / ideal : 0.0;
        fprintf(f, "%.10g %.10g %.10g %.10g %.10g %.10g\n", A->bins[b].r_center, pairs, sa, pd, ideal, gr);
    }

    fclose(f);
    return 0;
}
