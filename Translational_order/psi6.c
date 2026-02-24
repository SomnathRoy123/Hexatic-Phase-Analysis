/*
 * psi6.c
 *
 * Compute psi6 for each point given neighbor lists.
 * Uses minimum-image (mic_delta) for bond angles when use_pbc=true.
 */

#include "psi6.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

Complex *compute_psi6_from_neighbors(const Vec2Array *coms,
                                     const IntArray *neighbors,
                                     bool use_pbc,
                                     double box_x, double box_y)
{
    if(!coms || !neighbors){ 

        fprintf(stderr,"compute_psi6: invalid args\n");
        return NULL; 
    }
    int M = (int)coms->n;
    if(M <= 0) return NULL;

    Complex *psi = (Complex*)calloc((size_t)M, sizeof(Complex));
    if(!psi){ fprintf(stderr,"compute_psi6: OOM\n"); return NULL; }

    for(int i=0;i<M;i++){
        int nc = (int)neighbors[i].n;
        if(nc <= 0){ 
            psi[i].re = 0.0; 
            psi[i].im = 0.0; 
            continue; 
        }

        double sx = 0.0;
        double sy = 0.0;

        for(size_t k=0;k<neighbors[i].n;k++){
            int j = neighbors[i].data[k];
            if(j < 0 || j >= M) continue; /* defensive */

            double dx = coms->data[j].x - coms->data[i].x;
            double dy = coms->data[j].y - coms->data[i].y;
            if(use_pbc){
                dx = mic_delta(dx, box_x);
                dy = mic_delta(dy, box_y);
            }

            double theta = atan2(dy, dx);
            double ang6 = 6.0 * theta;
            sx += cos(ang6);
            sy += sin(ang6);
        }

        psi[i].re = sx / (double)nc;
        psi[i].im = sy / (double)nc;
    }

    return psi;
}


double compute_global_orientation_angle(const Complex *psi6, int M){
    if(!psi6 || M <= 0){
        fprintf(stderr, "compute_global_orientation_angle: invalid args\n");
        return 0.0;
    }

    double re_sum = 0.0;
    double im_sum = 0.0;
    for(int i=0;i<M;i++){
        re_sum += psi6[i].re;
        im_sum += psi6[i].im;
    }

    double invM = 1.0 / (double)M;
    double Theta6 = atan2(im_sum * invM, re_sum * invM);
    return Theta6 / 6.0;
}
