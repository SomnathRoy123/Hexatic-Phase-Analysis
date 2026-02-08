#include "utils.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h> /* for SIZE_MAX */

int v2a_init(Vec2Array *v){
    v->data = NULL;
    v->n = 0;
    v->cap = 0;
    return 0;
}

/* return 0 on success, -1 on OOM */
int v2a_push(Vec2Array *v, Vec2 p){
    if(v->n == v->cap){
        size_t newcap = v->cap ? v->cap * 2 : 128;
        if(newcap <= v->cap || newcap > (SIZE_MAX / sizeof(Vec2))) {
            /* overflow detected */
            return -1;
        }
        Vec2 *tmp = (Vec2*)realloc(v->data, newcap * sizeof(Vec2));
        if(!tmp){
            return -1; /* OOM */
        }
        v->data = tmp;
        v->cap = newcap;
    }
    v->data[v->n++] = p;
    return 0;
}

void v2a_free(Vec2Array *v){
    free(v->data);
    v->data = NULL;
    v->n = v->cap = 0;
}

/* IntArray */
int ia_init(IntArray *v){
    v->data = NULL;
    v->n = 0;
    v->cap = 0;
    return 0;
}

int ia_push(IntArray *v, int x){
    if(v->n == v->cap){
        size_t newcap = v->cap ? v->cap * 2 : 16;
        if(newcap <= v->cap || newcap > (SIZE_MAX / sizeof(int))){
            return -1;
        }
        int *tmp = (int*)realloc(v->data, newcap * sizeof(int));
        if(!tmp) return -1;
        v->data = tmp;
        v->cap = newcap;
    }
    v->data[v->n++] = x;
    return 0;
}

void ia_free(IntArray *v){
    free(v->data);
    v->data = NULL;
    v->n = v->cap = 0;
}

/* PBC helpers */
double mic_delta(double d, double L){
    if(L <= 0.0) return d;
    /* round is fine; */
    return d - L * round(d / L);
}

double wrap_pos(double x, double L){
    if(L <= 0.0) return x;
    double y = fmod(x, L);
    if(y < 0) y += L;
    return y;
}
