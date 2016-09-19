#ifndef __INTERVAL_TREE_H__
#define __INTERVAL_TREE_H__

#include <array>
#include <algorithm>

#include "global.h"
#include "itree_result.hpp"
#include "bound.h"


#define LAST_STOP(tree)   tree->by_stop.back()
#define FIRST_START(tree) tree->by_start.front()

#define T_SIZE(tree) tree->by_start.size()

#define T_START(tree, i) tree->by_start[i]
#define T_START_START(tree, i) tree->by_start[i]->pos[0]

#define T_STOP(tree, i) tree->by_stop[i]
#define T_STOP_STOP(tree, i) tree->by_stop[i]->pos[1]

#define R_SIZE(result) result->iv.size()

#define LEFT(tree) tree->children[0]
#define RIGHT(tree) tree->children[1]

typedef enum orientation {
    O_LEFT  = -1,
    O_ROOT  =  0,
    O_RIGHT =  1,
    O_UNSET =  2
} Orientation;

template <class T>
class IntervalTree
{
private:
    long get_center(std::vector<Bound*> intervals);

    /** Count the number of query-side intervals that overlap a search interval
     *
     * @param interval an interval on the IntervalTree
     * @param tree an IntervalTree
     */
    long count_overlaps(
        Bound* interval,
        IntervalTree* tree,
        long count
    );

    /** Count the number of query-side intervals that overlap a point
     *
     * @param point on the IntervalTree
     * @param tree an IntervalTree
     */
    long count_overlaps(
        long point,
        IntervalTree* tree,
        long count
    );

    /** Retrieve all target intervals overlapping a query interval
     *
     * First find all intervals query-side that overlap the query search interval.
     * Add all of these to the return IntervalResult object. If no intervals are
     * found, flag the search as 'inbetween' and return the first adjacent (by
     * start on target) block. Downstream functions can use this block as an anchor
     * and retrieve whatever else they need.
     *
     * @param query the query search interval
     * @param tree an IntervalTree
     */
    IntervalResult<T>* get_overlaps(
        Bound* interval,
        IntervalTree* tree,
        IntervalResult<T>* ir
    );

    /** Count the number of query-side intervals that overlap a point
     *
     * If the pnt if greater than the length of the sequence, then 0 is returned.
     *
     * @param pnt a point position on the IntervalTree
     * @param tree an IntervalTree
     */
    void get_overlaps(
        long point,
        IntervalTree* tree,
        IntervalResult<T>* ir
    );

    /**
     * Given a node containing one flank of a non-overlapping interval, find the other flank
     *                              ...
     *                             /
     *              _____________ B ______________
     *             /     |     |                  \
     *            /      |     | ----  i3           ...
     *           O       |     |  ---- i4
     *      ... / \      |     |   --- i5
     *             O     |     |
     *           /  \    |     |      Where:
     *       ...     A   |     |        - { i1, i2 } are intervals in A
     *           i1 ---- | --- |        - { i3, i4, i5 } are intervals in B
     *           i2 ---   query       Given A, adds i3 to result
     *
     * If starting from A, move up through parents until you reach a LEFT node.
     * Then ascend to the parent. This node is the lowest node on the right of A,
     * thus it will contain the nearest intervals.
     *
     */
    void set_nearest_opposing_interval(
        IntervalTree* tree,
        IntervalResult<T>* result,
        Pos position
    );

    // printer dispatcher
    void print(IntervalTree* n, int depth, char pos, int verbosity);
    void print_verbosity_1(IntervalTree* n, int depth, char pos);
    void print_verbosity_2(IntervalTree* n, int depth, char pos);
    void print_verbosity_3(IntervalTree* n, int depth, char pos);

protected:
    // the center position for this node
    long center;
    // all intervals that overlap the center, sorted by start position
    std::vector<T*> by_start;
    // all intervals that overlap the center, sorted by stop position
    std::vector<T*> by_stop;
    // Child nodes
    IntervalTree* children[2];
    // Parent
    IntervalTree* parent;
    // position relative to parent
    Orientation orientation;

public:
    IntervalTree(
        std::vector<T*> intervals,
        IntervalTree<T>* parent = NULL,
        Orientation orientation = O_ROOT
    );

    ~IntervalTree();

    void print(int verbosity);

    long count_overlaps(Bound* interval);

    long count_overlaps(long point);

    IntervalResult<T>* get_overlaps(Bound* interval);

    IntervalResult<T>* get_overlaps(long point);
};

#endif
