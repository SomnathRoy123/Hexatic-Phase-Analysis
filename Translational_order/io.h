#ifndef IO_H
#define IO_H

#include "utils.h"

/*
 * read_snapshot_xy
 *
 * Reads an ASCII snapshot file containing at least two columns: x y
 * Extra columns are ignored.
 * Lines starting with '#' or blank lines are skipped.
 *
 * On success:
 *   - fills Vec2Array *pos (must be uninitialized; function will call v2a_init)
 *   - returns 1
 *
 * On failure:
 *   - returns 0 (pos is left in a safe state)
 *
 * Caller must call v2a_free(pos) after use.
 */
int read_snapshot_xy(const char *path, Vec2Array *pos);

/*
 * extract_time_index
 *
 * Matches "time_<number>.dat" at the end of a path.
 * Returns integer index, or -1 if not matched.
 */
// int extract_time_index(const char *path);

/*
 * Natural sort comparison for paths containing time_<num>.dat
 * Use in qsort: qsort(paths, n, sizeof(char*), cmp_paths_by_time)
 */
// int cmp_paths_by_time(const void *a, const void *b);

#endif /* IO_H */
