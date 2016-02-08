#ifndef __INTERVAL_TREE_H__
#define __INTERVAL_TREE_H__

#include "ia.h"

struct IntervalTree {
    unsigned int center;
    IA * by_start;
    IA * by_stop;
    struct IntervalTree * l_child;
    struct IntervalTree * r_child;
};

struct IntervalTree * init_interval_tree();

void free_interval_tree(struct IntervalTree *);

#endif
