#ifndef UTILS_H
#define UTILS_H

#include <stddef.h>

/* --------- Basic 2D vector type --------- */
typedef struct { double x, y; } Vec2;

/* --------- Vec2 dynamic array ------------ */
typedef struct {
    Vec2 *data;
    size_t n, cap;
} Vec2Array;

int v2a_init(Vec2Array *v);
int v2a_push(Vec2Array *v, Vec2 p);
void v2a_free(Vec2Array *v);

/* --------- Int dynamic array ------------- */
typedef struct {
    int *data;
    size_t n, cap;
} IntArray;

int ia_init(IntArray *v);
int ia_push(IntArray *v, int x);
void ia_free(IntArray *v);

/* --------- PBC helpers ------------------- */
double mic_delta(double d, double L);
double wrap_pos(double x, double L);

/* --------- Misc helpers ------------------ */
static inline double sq(double x){ return x*x; }

#endif
