#ifndef __ITREE_H__
#define __ITREE_H__

#include "ia.h"
#include "interval.h"

#ifdef uint
#define uint unsigned int
#endif

/** An interval tree data structure allowing log(n) searches for overlapping
 * intervals.  */
struct IntervalTree {
    unsigned int center;
    IA * by_start;
    IA * by_stop;
    struct IntervalTree * l_child;
    struct IntervalTree * r_child;
};

struct IntervalTree * init_interval_tree();

void free_interval_tree(struct IntervalTree *);

void print_interval_tree(struct IntervalTree*, int verbosity);

struct IntervalTree * build_tree(IA*);

#endif
