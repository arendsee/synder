#ifndef __ITREE_NODE_H__
#define __ITREE_NODE_H__

#include <algorithm>

#include "global.h"
#include "interval.h"

#define LAST_STOP(tree)   tree->by_stop.back()
#define FIRST_START(tree) tree->by_start.front()

#define T_SIZE(tree) tree->by_start.size()

#define T_START(tree, i) tree->by_start[i]
#define T_START_START(tree, i) tree->by_start[i]->start

#define T_STOP(tree, i) tree->by_stop[i]
#define T_STOP_STOP(tree, i) tree->by_stop[i]->stop

#define R_SIZE(result) result->iv.size()

#define LEFT(tree) tree->children[0]
#define RIGHT(tree) tree->children[1]

typedef enum orientation {
    O_LEFT  = -1,
    O_ROOT  =  0,
    O_RIGHT =  1,
    O_UNSET =  2
} Orientation;

/** An interval tree class for log(n) searches for overlapping intervals */
class IntervalTreeNode {
    private:
        long get_center(std::vector<Interval*> intervals);

    public:
        // the center position for this node
        long center;
        // all intervals that overlap the center, sorted by start position
        std::vector<Interval*> by_start;
        // all intervals that overlap the center, sorted by stop position
        std::vector<Interval*> by_stop;
        // Child nodes
        IntervalTreeNode * children[2];
        // Parent
        IntervalTreeNode * parent;
        // position relative to parent
        Orientation orientation;

        /** Recursive constructor */
        IntervalTreeNode(
            std::vector<Interval*> intervals,
            IntervalTreeNode * parent = NULL,
            Orientation orientation = O_ROOT
        );

        ~IntervalTreeNode();

        // Wrapper for the treal recursive printer
        void print(int verbosity);
        void print(IntervalTreeNode * n, int depth, char pos, int verbosity);
        void print_verbosity_1(IntervalTreeNode * n, int depth, char pos);
        void print_verbosity_2(IntervalTreeNode * n, int depth, char pos);
        void print_verbosity_3(IntervalTreeNode * n, int depth, char pos);
};

#endif
