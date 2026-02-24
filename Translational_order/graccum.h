#ifndef GRACCUM_H
#define GRACCUM_H

#include "utils.h"
#include <stdbool.h>

typedef struct GrAccum GrAccum;

GrAccum *graccum_create(double dr);
void graccum_free(GrAccum *A);

void graccum_accumulate(GrAccum *A,
                        const Vec2Array *coms,
                        bool use_pbc,
                        double box_x,
                        double box_y);

/* Estimate lattice spacing a from the first peak in averaged g(r).
 * Returns >0 on success, <=0 when no valid peak found.
 */
double graccum_estimate_first_peak_a(const GrAccum *A);

int graccum_write(GrAccum *A,
                  const char *outpath,
                  int t0,
                  int t1,
                  bool use_pbc,
                  double box_x,
                  double box_y);

#endif
