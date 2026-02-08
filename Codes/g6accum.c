/*
 * g6accum.c
 *
 * Accumulator for hexatic correlation g6(r).
 *
 * Usage:
 *   G6Accum *A = g6accum_create(dr);
 *   g6accum_accumulate(A, &coms, psi6, USE_PBC, BOX_X, BOX_Y);  // per snapshot
 *   g6accum_write(A, "g6_avg_time_100_200.dat", 100,200, LBOND, USE_PBC, BOX_X, BOX_Y);
 *   g6accum_free(A);
 */

#include "g6accum.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>



typedef struct {
    double r_center;
    double re_sum;
    double im_sum;
    long   pair_count;
} G6Bin;

struct G6Accum {
    G6Bin *bins;
    int    nbins;
    double dr;
};

/* Create accumulator */
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

/* Free accumulator */
void g6accum_free(G6Accum *A){
    if(!A) return;
    free(A->bins);
    A->bins = NULL;
    A->nbins = 0;
    free(A);
}

/* Ensure we have bins up to index bmax (inclusive) */
static void g6accum_ensure_bins(G6Accum *A, int bmax){
    if(bmax < 0) return;
    if(bmax < A->nbins) return;
    int new_n = bmax + 1;
    G6Bin *nb = (G6Bin*)realloc(A->bins, (size_t)new_n * sizeof(G6Bin));
    if(!nb){ 
        fprintf(stderr,"g6accum_ensure_bins: OOM\n"); 
        exit(1); 
    }
    /* initialize newly allocated bins */
    for(int b = A->nbins; b < new_n; ++b){
        nb[b].r_center = (b + 0.5) * A->dr;
        nb[b].re_sum = 0.0;
        nb[b].im_sum = 0.0;
        nb[b].pair_count = 0;
    }
    A->bins = nb;
    A->nbins = new_n;
}

/* Accumulate contributions from a snapshot.
 * For each unordered pair i<j of COMs:
 *   r = distance(coms[i], coms[j]) (use mic_delta for PBC),
 *   bin = floor(r / dr),
 *   add Re(psi_i * conj(psi_j)) and Im(...)
 *   increment pair_count
 */
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

    /* First pass: find rmax to know how many bins we need */
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
            if(r2 > rmax*rmax) rmax = sqrt(r2);
        }
    }
    int bmax = (int)floor(rmax / A->dr);
    g6accum_ensure_bins(A, bmax);

    /* Accumulate pair contributions */
    for(int i=0;i<M-1;i++){
        for(int j=i+1;j<M;j++){
            double dx = coms->data[j].x - coms->data[i].x;
            double dy = coms->data[j].y - coms->data[i].y;
            if(use_pbc){
                dx = mic_delta(dx, box_x);
                dy = mic_delta(dy, box_y);
            }
            double r = sqrt(dx*dx + dy*dy);
            int b = (int)floor(r / A->dr);
            if(b < 0 || b >= A->nbins) continue;

            /* (a + i b)*(c - i d) = (ac + bd) + i (bc - ad) */
            double a = psi6[i].re;
            double b_im = psi6[i].im;
            double c = psi6[j].re;
            double d = psi6[j].im;

            A->bins[b].re_sum += a * c + b_im * d;
            A->bins[b].im_sum += b_im * c - a * d;
            A->bins[b].pair_count += 1;
        }
    }
}

/* Write averaged g6(r) file. Format:
 * r_center  Re[g6(r)]  Im[g6(r)]  |g6(r)|  pair_count
 */
int g6accum_write(G6Accum *A,
                  const char *outpath,
                  int t0, int t1,
                  double lbond,
                  bool use_pbc,
                  double box_x,
                  double box_y)
{
    if(!A || !outpath){ fprintf(stderr,"g6accum_write: invalid args\n"); return 1; }
    FILE *f = fopen(outpath, "w");
    if(!f){ fprintf(stderr,"g6accum_write: cannot open %s: %s\n", outpath, strerror(errno)); return 2; }

    fprintf(f, "# Averaged g6(r) over snapshots time_%d .. time_%d\n", t0, t1);
    fprintf(f, "# Columns: r_center  Re[g6(r)]  Im[g6(r)]  |g6(r)|  pair_count\n");
    fprintf(f, "# Params: dr = %.8g  lbond = %.8g  USE_PBC = %s\n", A->dr, lbond, use_pbc ? "true" : "false");
    if(use_pbc){
        fprintf(f, "# Box dims: %.8g x %.8g\n", box_x, box_y);
    }

    for(int b=0;b<A->nbins;b++){
        long cnt = A->bins[b].pair_count;
        if(cnt <= 0) continue;
        double re = A->bins[b].re_sum / (double)cnt;
        double im = A->bins[b].im_sum / (double)cnt;
        double mag = sqrt(re*re + im*im);
        fprintf(f, "%.8f %.10e %.10e %.10e %ld\n", A->bins[b].r_center, re, im, mag, cnt);
    }

    fclose(f);
    return 0;
}
