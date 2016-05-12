#ifndef __ITREE_H__
#define __ITREE_H__

#include "ia.h"
#include "interval.h"

#ifdef uint
#define uint unsigned int
#endif

/** An interval tree data structure for log(n) searches for overlapping intervals */
struct IntervalTree {
    // the center position for this node
    unsigned int center;
    // all intervals that overlap the center, sorted by start position
    IA * by_start;
    // all intervals that overlap the center, sorted by stop position
    IA * by_stop;
    // Child node with center less than this center
    struct IntervalTree * l_child;
    // Child node with center greater than this center
    struct IntervalTree * r_child;
};

struct IntervalTree * init_interval_tree();

void free_interval_tree(struct IntervalTree *);

void print_interval_tree(struct IntervalTree*, int verbosity);

struct IntervalTree * build_tree(IA*);

#endif
