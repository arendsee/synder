#ifndef __ITREE_H__
#define __ITREE_H__

#include <array>

#include "global.h"
#include "interval.h"
#include "itree_node.h"
#include "itree_result.h"

class IntervalTree
{
private:
    size_t size;
    Interval * interval_pool;
    IntervalTreeNode * root;

    /** Count the number of query-side intervals that overlap a search interval
     *
     * @param interval an interval on the IntervalTree
     * @param tree an IntervalTree
     */
    long count_overlaps(
        Interval * interval,
        IntervalTreeNode * tree,
        long count
    );


    /** Count the number of query-side intervals that overlap a point
     *
     * @param point on the IntervalTree
     * @param tree an IntervalTree
     */
    long count_overlaps(
        long point,
        IntervalTreeNode * tree,
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
    IntervalResult * get_overlaps(
        Interval * interval,
        IntervalTreeNode * tree,
        IntervalResult * ir
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
        IntervalTreeNode * tree,
        IntervalResult * ir
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
        IntervalTreeNode * tree,
        IntervalResult * result,
        Pos position
    );

public:
    IntervalTree(Interval * ipool, size_t ipool_size);

    ~IntervalTree();

    void print(int verbosity);

    /** Wrapper for search.h::count_overlaps */
    long count_overlaps(Interval * interval);

    /** Wrapper for search.h::count_overlaps */
    long count_overlaps(long point);

    /** Wrapper for search.h::get_overlaps */
    IntervalResult * get_overlaps(Interval * interval);

    /** Wrapper for search.h::get_overlaps */
    IntervalResult * get_overlaps(long point);
};

#endif
