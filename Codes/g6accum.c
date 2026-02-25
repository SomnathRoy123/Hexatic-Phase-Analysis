/*
 * g6accum.c
 *
 * Accumulator for hexatic correlation g6(r).
 */

#include "g6accum.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>

typedef struct {
    double r_center;
    double re_sum;        /* sum over per-snapshot bin means */
    double im_sum;        /* sum over per-snapshot bin means */
    long   pair_count;    /* total contributing pairs (diagnostic) */
    long   sample_count;  /* number of snapshots contributing to this bin */
} G6Bin;

struct G6Accum {
    G6Bin *bins;
    int    nbins;
    double dr;
};

G6Accum *g6accum_create(double dr){
    if(dr <= 0.0){
        fprintf(stderr, "g6accum_create: dr must be > 0\n");
        return NULL;
    }
    G6Accum *A = (G6Accum*)malloc(sizeof(G6Accum));
    if(!A){
        fprintf(stderr,"g6accum_create: OOM\n");
        return NULL;
    }
    A->bins = NULL;
    A->nbins = 0;
    A->dr = dr;
    return A;
}

void g6accum_free(G6Accum *A){
    if(!A) return;
    free(A->bins);
    A->bins = NULL;
    A->nbins = 0;
    free(A);
}

static void g6accum_ensure_bins(G6Accum *A, int bmax){
    if(bmax < 0) return;
    if(bmax < A->nbins) return;
    int new_n = bmax + 1;
    G6Bin *nb = (G6Bin*)realloc(A->bins, (size_t)new_n * sizeof(G6Bin));
    if(!nb){
        fprintf(stderr,"g6accum_ensure_bins: OOM\n");
        exit(1);
    }
    for(int b = A->nbins; b < new_n; ++b){
        nb[b].r_center = (b + 0.5) * A->dr;
        nb[b].re_sum = 0.0;
        nb[b].im_sum = 0.0;
        nb[b].pair_count = 0;
        nb[b].sample_count = 0;
    }
    A->bins = nb;
    A->nbins = new_n;
}

void g6accum_accumulate(G6Accum *A,
                        const Vec2Array *coms,
                        const Complex *psi6,
                        bool use_pbc,
                        double box_x,
                        double box_y)
{
    if(!A || !coms || !psi6) return;
    int M = (int)coms->n;
    if(M < 2) return;

    /*
     * In periodic boxes, isotropic radial correlations are only reliable up to
     * half the smallest box length. Larger r values can create artificial tails.
     */
    double max_valid_r = INFINITY;
    if(use_pbc){
        max_valid_r = 0.5 * fmin(box_x, box_y);
        if(max_valid_r <= 0.0) return;
    }

    double rmax = 0.0;
    for(int i=0;i<M-1;i++){
        for(int j=i+1;j<M;j++){
            double dx = coms->data[j].x - coms->data[i].x;
            double dy = coms->data[j].y - coms->data[i].y;
            if(use_pbc){
                dx = mic_delta(dx, box_x);
                dy = mic_delta(dy, box_y);
            }
            double r2 = dx*dx + dy*dy;
            if(use_pbc && r2 > max_valid_r*max_valid_r) continue;
            if(r2 > rmax*rmax) rmax = sqrt(r2);
        }
    }

    int bmax = (int)floor(rmax / A->dr);
    g6accum_ensure_bins(A, bmax);

    /*
     * Per-snapshot bin means (frame-equal averaging): each snapshot contributes
     * one average per bin, preventing overweighting by snapshots with more pairs.
     */
    double *snap_re = (double*)calloc((size_t)A->nbins, sizeof(double));
    double *snap_im = (double*)calloc((size_t)A->nbins, sizeof(double));
    long   *snap_cnt = (long*)calloc((size_t)A->nbins, sizeof(long));
    if(!snap_re || !snap_im || !snap_cnt){
        free(snap_re);
        free(snap_im);
        free(snap_cnt);
        return;
    }

    for(int i=0;i<M-1;i++){
        for(int j=i+1;j<M;j++){
            double dx = coms->data[j].x - coms->data[i].x;
            double dy = coms->data[j].y - coms->data[i].y;
            if(use_pbc){
                dx = mic_delta(dx, box_x);
                dy = mic_delta(dy, box_y);
            }
            double r = sqrt(dx*dx + dy*dy);
            if(use_pbc && r > max_valid_r) continue;
            int b = (int)floor(r / A->dr);
            if(b < 0 || b >= A->nbins) continue;

            /* (a + i b)*(c - i d) = (ac + bd) + i (bc - ad) */
            double a = psi6[i].re;
            double b_im = psi6[i].im;
            double c = psi6[j].re;
            double d = psi6[j].im;

            snap_re[b] += a * c + b_im * d;
            snap_im[b] += b_im * c - a * d;
            snap_cnt[b] += 1;
        }
    }

    for(int b=0;b<A->nbins;b++){
        if(snap_cnt[b] <= 0) continue;
        A->bins[b].re_sum += snap_re[b] / (double)snap_cnt[b];
        A->bins[b].im_sum += snap_im[b] / (double)snap_cnt[b];
        A->bins[b].pair_count += snap_cnt[b];
        A->bins[b].sample_count += 1;
    }

    free(snap_re);
    free(snap_im);
    free(snap_cnt);
}

int g6accum_write(G6Accum *A,
                  const char *outpath,
                  int t0, int t1,
                  double lbond,
                  bool use_pbc,
                  double box_x,
                  double box_y)
{
    if(!A || !outpath){
        fprintf(stderr,"g6accum_write: invalid args\n");
        return 1;
    }
    FILE *f = fopen(outpath, "w");
    if(!f){
        fprintf(stderr,"g6accum_write: cannot open %s: %s\n", outpath, strerror(errno));
        return 2;
    }

    fprintf(f, "# Averaged g6(r) over snapshots time_%d .. time_%d\n", t0, t1);
    fprintf(f, "# Columns: r_center  Re[g6(r)]  Im[g6(r)]  |g6(r)|  n_snapshots  pair_count_total\n");
    fprintf(f, "# Params: dr = %.8g  lbond = %.8g  USE_PBC = %s\n", A->dr, lbond, use_pbc ? "true" : "false");
    if(use_pbc){
        double rcut = 0.5 * fmin(box_x, box_y);
        fprintf(f, "# Box dims: %.8g x %.8g\n", box_x, box_y);
        fprintf(f, "# PBC radial cutoff applied: r <= %.8g\n", rcut);
    }

    for(int b=0;b<A->nbins;b++){
        long ns = A->bins[b].sample_count;
        long cnt = A->bins[b].pair_count;
        if(ns <= 0) continue;
        double re = A->bins[b].re_sum / (double)ns;
        double im = A->bins[b].im_sum / (double)ns;
        double mag = sqrt(re*re + im*im);
        fprintf(f, "%.8f %.10e %.10e %.10e %ld %ld\n",
                A->bins[b].r_center, re, im, mag, ns, cnt);
    }

    fclose(f);
    return 0;
}
