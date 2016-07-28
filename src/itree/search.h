#ifndef __SEARCH_H__
#define __SEARCH_H__

#include <stdbool.h>

#include "ia.h"
#include "iv.h"
#include "interval.h"
#include "itree.h"

#ifndef uint
#define uint unsigned int
#endif

#define LAST_STOP(tree)   tree->by_stop->v[tree->by_stop->size-1]
#define FIRST_START(tree) tree->by_start->v[0]

#define T_SIZE(tree) tree->by_start->size

#define T_START(tree, i) tree->by_start->v[i]
#define T_START_START(tree, i) tree->by_start->v[i].start

#define T_STOP(tree, i) tree->by_stop->v[i]
#define T_STOP_STOP(tree, i) tree->by_stop->v[i].stop

#define R_SIZE(result) result->iv->size

/** A container for results from searches for overlaps
 *
 * IV * iv        - a list of the blocks overlapping the query
 * bool inbetween - TRUE if the query overlaps nothing
 */
typedef struct {
    IV * iv;
    bool inbetween;
    bool leftmost;
    bool rightmost;
} IntervalResult;


/** Count the number of query-side intervals that overlap a point
 *
 * If the pnt if greater than the length of the sequence, then 0 is returned.
 *
 * @param pnt a point position on the IntervalTree
 * @param tree an IntervalTree
 */
uint count_point_overlaps(uint pnt, IntervalTree * tree);


/** Count the number of query-side intervals that overlap a search interval
 *
 * @param interval an interval on the IntervalTree
 * @param tree an IntervalTree
 */
uint count_interval_overlaps(Interval * interval, IntervalTree * tree);


/** Allocate memory and set defaults for an IntervalResults struct */
IntervalResult * init_IntervalResult();


/** Free memory of an IntervalResults struct */
void free_IntervalResult(IntervalResult *);

/** Print interval result */
void print_IntervalResult(IntervalResult *);


/** Retrieve all target intervals overlapping a query point
 *
 * See notes on get_point_overlaps
 */
IntervalResult * get_point_overlaps(uint, IntervalTree *);

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
IntervalResult * get_interval_overlaps(Interval *, IntervalTree *);

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
void set_nearest_opposing_interval(IntervalTree * tree, IntervalResult * result, Pos position);

#endif
