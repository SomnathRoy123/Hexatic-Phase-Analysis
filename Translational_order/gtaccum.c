#include "gtaccum.h"
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
    double ct_sum;
    long pair_count;
} GtBin;

struct GtAccum {
    GtBin *bins;
    int nbins;
    double dr;
};

GtAccum *gtaccum_create(double dr){
    if(dr <= 0.0){
        fprintf(stderr, "gtaccum_create: dr must be > 0\n");
        return NULL;
    }
    GtAccum *A = (GtAccum*)malloc(sizeof(GtAccum));
    if(!A){
        fprintf(stderr, "gtaccum_create: OOM\n");
        return NULL;
    }
    A->bins = NULL;
    A->nbins = 0;
    A->dr = dr;
    return A;
}

void gtaccum_free(GtAccum *A){
    if(!A) return;
    free(A->bins);
    free(A);
}

static void ensure_bins(GtAccum *A, int bmax){
    if(bmax < 0 || bmax < A->nbins) return;
    int new_n = bmax + 1;
    GtBin *nb = (GtBin*)realloc(A->bins, (size_t)new_n * sizeof(GtBin));
    if(!nb){
        fprintf(stderr, "gtaccum_ensure_bins: OOM\n");
        exit(1);
    }
    for(int b = A->nbins; b < new_n; ++b){
        nb[b].r_center = (b + 0.5) * A->dr;
        nb[b].ct_sum = 0.0;
        nb[b].pair_count = 0;
    }
    A->bins = nb;
    A->nbins = new_n;
}

void gtaccum_accumulate(GtAccum *A,
                        const Vec2Array *coms,
                        double theta_G,
                        double a_lattice,
                        bool use_pbc,
                        double box_x,
                        double box_y){
    if(!A || !coms) return;
    int M = (int)coms->n;
    if(M < 2) return;
    if(a_lattice <= 0.0){
        fprintf(stderr, "gtaccum_accumulate: a_lattice must be > 0\n");
        return;
    }
    if(use_pbc && (box_x <= 0.0 || box_y <= 0.0)){
        fprintf(stderr, "gtaccum_accumulate: invalid PBC box\n");
        return;
    }

    double G_mag = 4.0 * M_PI / (a_lattice * sqrt(3.0));
    double Gx[6], Gy[6];
    for(int n = 0; n < 6; ++n){
        double ang = theta_G + n * (M_PI / 3.0);
        Gx[n] = G_mag * cos(ang);
        Gy[n] = G_mag * sin(ang);
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
            if(b < 0 || b >= A->nbins) continue;

            double ct = 0.0;
            for(int n = 0; n < 6; ++n){
                ct += cos(Gx[n] * dx + Gy[n] * dy);
            }
            ct /= 6.0;

            A->bins[b].ct_sum += ct;
            A->bins[b].pair_count += 1;
        }
    }
}

int gtaccum_write(GtAccum *A,
                  const char *outpath,
                  int t0,
                  int t1,
                  double a_lattice,
                  bool use_pbc,
                  double box_x,
                  double box_y){
    if(!A || !outpath) return 1;

    FILE *f = fopen(outpath, "w");
    if(!f){
        fprintf(stderr, "gtaccum_write: cannot open %s: %s\n", outpath, strerror(errno));
        return 2;
    }

    fprintf(f, "# gT(r) average over time_%d..time_%d\n", t0, t1);
    fprintf(f, "# columns: r_center gT_avg pair_count\n");
    fprintf(f, "# dr=%.8g a_lattice=%.8g use_pbc=%s\n", A->dr, a_lattice, use_pbc ? "true" : "false");
    if(use_pbc) fprintf(f, "# box: %.8g %.8g\n", box_x, box_y);

    for(int b = 0; b < A->nbins; ++b){
        double gt = (A->bins[b].pair_count > 0) ? (A->bins[b].ct_sum / (double)A->bins[b].pair_count) : 0.0;
        fprintf(f, "%.10g %.10g %ld\n", A->bins[b].r_center, gt, A->bins[b].pair_count);
    }

    fclose(f);
    return 0;
}
