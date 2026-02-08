#ifndef G6ACCUM_H
#define G6ACCUM_H

#include "utils.h"   /* Vec2Array, mic_delta */
#include "psi6.h"    /* Complex */
#include <stdbool.h>

/* Opaque accumulator */
typedef struct G6Accum G6Accum;

/* Create/destroy */
G6Accum *g6accum_create(double dr);
void g6accum_free(G6Accum *A);

/* Accumulate one snapshot's contributions.
 * - coms: Vec2Array of M cluster COMs
 * - psi6: Complex array length M (psi6 at each COM)
 * - use_pbc / box dims: for distance calculation (if use_pbc is true, box_x/box_y must be >0)
 */
void g6accum_accumulate(G6Accum *A,
                        const Vec2Array *coms,
                        const Complex *psi6,
                        bool use_pbc,
                        double box_x,
                        double box_y);

/* Write averaged g6 file:
 * - outpath: path to output file to create (filename, not directory)
 * - t0,t1: time index range used in header
 * - lbond, use_pbc etc appear in header
 * Returns 0 on success, non-zero on failure.
 */
int g6accum_write(G6Accum *A,
                  const char *outpath,
                  int t0, int t1,
                  double lbond,
                  bool use_pbc,
                  double box_x,
                  double box_y);

#endif /* G6ACCUM_H */
