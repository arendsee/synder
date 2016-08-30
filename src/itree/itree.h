#ifndef __ITREE_H__
#define __ITREE_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

#include "interval.h"
#include "ia.h"
#include "iv.h"

typedef enum orientation {
    O_LEFT = -1,
    O_ROOT = 0,
    O_RIGHT = 1,
    O_UNSET = 99
} Orientation;

#define LEFT(tree) tree->children[0]
#define RIGHT(tree) tree->children[1]

/** An interval tree data structure for log(n) searches for overlapping intervals */
typedef struct IntervalTree {
    // the center position for this node
    long center;
    // all intervals that overlap the center, sorted by start position
    IA * by_start;
    // all intervals that overlap the center, sorted by stop position
    IA * by_stop;
    // Child nodes
    struct IntervalTree * children[2];
    // Parent
    struct IntervalTree * parent;
    // position relative to parent
    Orientation orientation;
} IntervalTree;

IntervalTree * init_IntervalTree();

void free_IntervalTree(IntervalTree *);

void print_IntervalTree(IntervalTree*, int verbosity);

IntervalTree * build_tree(IA*);

#endif
