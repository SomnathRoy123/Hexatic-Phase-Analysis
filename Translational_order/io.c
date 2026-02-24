#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* --------------------- read_snapshot_xy --------------------- */
int read_snapshot_xy(const char *path, Vec2Array *pos)
{
    if(!path || !pos){
        fprintf(stderr, "read_snapshot_xy: invalid arguments\n");
        return 0;
    }

    FILE *fp = fopen(path, "r");
    if(!fp){
        fprintf(stderr, "read_snapshot_xy: cannot open %s\n", path);
        return 0;
    }

    // v2a_init(pos);

    char line[4096];
    while(fgets(line, sizeof(line), fp)){
        /* skip comments and blank lines */
        if(line[0] == '#' || line[0] == '\n') continue;

        double x, y, z;
        if(sscanf(line, "%lf %lf  %lf", &x, &y , &z) == 3){
            v2a_push(pos, (Vec2){x, y});
        }
        /* else ignore malformed lines */
    }

    fclose(fp);
    return 1;
}

/* --------------------- extract_time_index --------------------- */
/* Extracts trailing time_<index>.dat */
int extract_time_index(const char *path)
{
    if(!path) return -1;

    /* find last '/' */
    const char *name = strrchr(path, '/');
    name = name ? name + 1 : path;

    int idx = -1;
    if(sscanf(name, "time_%d.dat", &idx) == 1) return idx;
    return -1;
}

/* --------------------- cmp_paths_by_time ---------------------- */
int cmp_paths_by_time(const void *a, const void *b)
{
    const char *pa = *(const char * const *)a;
    const char *pb = *(const char * const *)b;

    int ia = extract_time_index(pa);
    int ib = extract_time_index(pb);

    if(ia < ib) return -1;
    if(ia > ib) return 1;

    /* fallback: lexical */
    return strcmp(pa, pb);
}
