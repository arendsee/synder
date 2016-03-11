#ifndef __SEARCH_H__
#define __SEARCH_H__

#include <stdbool.h>

#include "ia.h"
#include "iv.h"
#include "interval.h"
#include "itree.h"

#ifdef uint
#define uint unsigned int
#endif

#define LAST_STOP(tree)   tree->by_stop->v[tree->by_stop->size-1]
#define FIRST_START(tree) tree->by_start->v[0]

#define T_SIZE(tree) tree->by_start->size

#define T_START(tree, i) tree->by_start->v[i]
#define T_START_START(tree, i) tree->by_start->v[i].start

#define T_STOP(tree, i) tree->by_stop->v[i]
#define T_STOP_STOP(tree, i) tree->by_stop->v[i].stop

#define RIGHT(tree) tree->r_child
#define LEFT(tree) tree->l_child

#define R_SIZE(result) result->iv->size

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
