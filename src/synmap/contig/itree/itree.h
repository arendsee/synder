#ifndef __ITREE_H__
#define __ITREE_H__

#include "ia.h"
#include "interval.h"
#include "interval-tree.h"

#ifdef uint
#define uint unsigned int
#endif

Pos point_overlap(uint, Interval);

struct IntervalTree * build_tree(IA*);

void print_interval_tree(struct IntervalTree*, int verbosity);

#endif
