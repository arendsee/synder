#ifndef __SEARCH_H__
#define __SEARCH_H__

#include <stdbool.h>

#include "ia.h"
#include "iv.h"
#include "interval.h"
#include "interval-tree.h"

#ifdef uint
#define uint unsigned int
#endif

#define LAST_STOP(tree)   tree->by_stop->v[tree->by_stop->size]
#define FIRST_START(tree) tree->by_start->v[0]

typedef struct {
    IV * iv;
    bool inbetween;
} IntervalResult;

uint count_point_overlaps(uint, struct IntervalTree *);

uint count_interval_overlaps(Interval *, struct IntervalTree *);

IntervalResult * init_IntervalResult(size_t initial_iv_size);

void free_IntervalResult(IntervalResult *);

IntervalResult * get_point_overlaps(uint, struct IntervalTree *);

IntervalResult * get_interval_overlaps(Interval *, struct IntervalTree *);

#endif
