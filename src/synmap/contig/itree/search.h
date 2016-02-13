#ifndef __SEARCH_H__
#define __SEARCH_H__

#include "ia.h"
#include "interval.h"
#include "interval-tree.h"

#ifdef uint
#define uint unsigned int
#endif

uint count_point_overlaps(uint, struct IntervalTree *);

uint count_interval_overlaps(Interval *, struct IntervalTree *);

IA * get_point_overlaps(uint, struct IntervalTree *);

IA * get_interval_overlaps(Interval *, struct IntervalTree *);

#endif
