#ifndef __INTERVAL_TREE_H__
#define __INTERVAL_TREE_H__

#include <array>
#include <algorithm>

#include "global.h"
#include "interval_result.hpp"


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

    static bool cmp_start (T* a, T* b) { return ( a->pos[0] < b->pos[0] ); }
    static bool cmp_stop  (T* a, T* b) { return ( a->pos[1] < b->pos[1] ); }

    /*
     * Select a point at the center of the middle interval.
     * This guarantees at least one interval overlaps each node.
     * If the intervals are sorted, it also favors (but doesn't guarantee) a
     * balanced tree.
     */
    long get_center(std::vector<T*> v)
    {
        // get the central index
        long i = v.size() / 2;
        // get the center point on this index
        long x = (STOP(v[i]) - START(v[i])) / 2 + START(v[i]);
        return x;
    }

    /** Count the number of query-side intervals that overlap a search interval
     *
     * @param interval an interval on the IntervalTree
     * @param tree an IntervalTree
     */
    template <class U>
    long count_overlaps(U* inv, IntervalTree<T>* tree, long count)
    {
        if (tree == nullptr)
            return count;
        switch (inv->position_relative_to(tree->center)) {
            case lo:
                for (long long i = T_SIZE(tree) - 1; i >= 0 ; i--) {
                    if (START(inv) <= T_STOP_STOP(tree, i)) {
                        count++;
                    } else {
                        break;
                    }
                }
                return count_overlaps(inv, RIGHT(tree), count);
            case hi:
                for (size_t i = 0; i < T_SIZE(tree); i++) {
                    if (STOP(inv) >= T_START_START(tree, i)) {
                        count++;
                    } else {
                        break;
                    }
                }
                return count_overlaps(inv, LEFT(tree), count);
            default:
                count += T_SIZE(tree);
                return count_overlaps(inv, RIGHT(tree),
                                      count_overlaps(inv, LEFT(tree), count));
        }
    }

    /** Count the number of query-side intervals that overlap a point
     *
     * @param point on the IntervalTree
     * @param tree an IntervalTree
     */
    long count_overlaps(long pnt, IntervalTree<T>* tree, long count)
    {
        if (pnt >= tree->center) {
            for (long long i = T_SIZE(tree) - 1; i >= 0 ; i--) {
                if (pnt <= T_STOP_STOP(tree, i)) {
                    count++;
                } else {
                    break;
                }
            }
            if (RIGHT(tree) != nullptr)
                return count_overlaps(pnt, RIGHT(tree), count);
        } else {
            for (size_t i = 0; i < T_SIZE(tree); i++) {
                if (pnt >= T_START_START(tree, i)) {
                    count++;
                } else {
                    break;
                }
            }
            if (LEFT(tree) != nullptr)
                return count_overlaps(pnt, LEFT(tree), count);
        }
        return count;
    }

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
    template <class U>
    IntervalResult<T>* get_overlaps(U* inv, IntervalTree<T>* tree, IntervalResult<T>* results)
    {
        if (tree == nullptr)
            return results;
        switch (inv->position_relative_to(tree->center)) {
            case lo: // center lower than interval start
                for (long long i = T_SIZE(tree) - 1; i >= 0 ; i--) {
                    if (START(inv) <= T_STOP_STOP(tree, i)) {
                        results->iv.push_back(T_STOP(tree, i));
                    } else {
                        break;
                    }
                }
                // Reach a leaf and still have no overlaps
                if (RIGHT(tree) == nullptr && results->iv.size() == 0) {
                    results->iv.push_back(LAST_STOP(tree));
                    // get nearest interval on the other side
                    set_nearest_opposing_interval(tree, results, hi);
                }
                return get_overlaps(inv, RIGHT(tree), results);
            case hi:
                for (size_t i = 0; i < T_SIZE(tree); i++) {
                    if (STOP(inv) >= T_START_START(tree, i)) {
                        results->iv.push_back(tree->by_start[i]);
                    } else {
                        break;
                    }
                }
                // Reach a leaf and still have no overlaps
                if (LEFT(tree) == nullptr && results->iv.size() == 0) {
                    results->iv.push_back(FIRST_START(tree));
                    // get nearest interval on the other side
                    set_nearest_opposing_interval(tree, results, lo);
                }
                return get_overlaps(inv, LEFT(tree), results);
            default: // in
                results->iv.insert(
                    results->iv.end(),
                    tree->by_start.begin(),
                    tree->by_start.end()
                );
                // join_IV(results->iv, tree->by_start);
                return get_overlaps(inv, RIGHT(tree),
                                    get_overlaps(inv, LEFT(tree), results));
        };
    }

    /** Count the number of query-side intervals that overlap a point
     *
     * If the pnt if greater than the length of the sequence, then 0 is returned.
     *
     * @param pnt a point position on the IntervalTree
     * @param tree an IntervalTree
     */
    void get_point_overlaps(long pnt, IntervalTree<T>* tree, IntervalResult<T>* results)
    {
        if (pnt >= tree->center) {
            for (long long i = T_SIZE(tree) - 1; i >= 0 ; i--) {
                if (pnt <= T_STOP_STOP(tree, i)) {
                    results->iv.push_back(T_STOP(tree, i));
                } else {
                    break;
                }
            }
            if (RIGHT(tree) != nullptr) {
                get_overlaps(pnt, RIGHT(tree), results);
            } else if (R_SIZE(results) != 0) {
                results->iv.push_back(LAST_STOP(tree));
                set_nearest_opposing_interval(tree, results, hi);
            }
        } else {
            for (size_t i = 0; i < T_SIZE(tree); i++) {
                if (pnt >= T_START_START(tree, i)) {
                    results->iv.push_back(T_START(tree, i));
                } else {
                    break;
                }
            }
            if (LEFT(tree) != nullptr) {
                get_overlaps(pnt, LEFT(tree), results);
            } else if (R_SIZE(results) != 0) {
                results->iv.push_back(FIRST_START(tree));
                set_nearest_opposing_interval(tree, results, lo);
            }
        }
    }

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
    void set_nearest_opposing_interval(IntervalTree<T>* node, IntervalResult<T>* results, Pos pos)
    {

        Orientation l_orientation = node->orientation;

        if ((pos == hi && l_orientation == O_LEFT) ||
                (pos == lo && l_orientation == O_RIGHT)) {
            node = node->parent;
        } else {
            while (node->orientation != O_ROOT && node->orientation == l_orientation) {
                node = node->parent;
            }
            node = node->parent;
        }

        if (node == nullptr) {
            if (pos == lo) {
                results->leftmost = true;
            } else {
                results->rightmost = true;
            }
        } else {
            results->inbetween = true;
            if (l_orientation == O_RIGHT) {
                results->iv.push_back(node->by_start.front());
            } else if (l_orientation == O_LEFT) {
                results->iv.push_back(node->by_stop[node->by_stop.size() - 1]);
            }
        }
    }

    // printer dispatcher
    /* write tree and center */
    void print_verbosity_1(IntervalTree<T>* n, int depth, char pos)
    {
        printf("%*d - %c%zu\n", depth * 2, depth, pos, n->center);
    }

    /* write tree, center, and start-sorted */
    void print_verbosity_2(IntervalTree<T>* n, int depth, char pos)
    {
        printf("%*d   %*s\t%c%zu:",
               depth * 2, depth,
               10 - depth * 2, "", pos, n->center);
        for(size_t i = 0; i < n->by_start.size(); i++) {
            printf("(%zu,%zu) ",
                   START(n->by_start[i]),
                   STOP(n->by_start[i]));
        }
        printf("\n");
    }


    /* write start- and stop-sorted vectors for each node */
    void print_verbosity_3(IntervalTree<T>* n, int depth, char pos)
    {
        print_verbosity_1(n, depth, pos);
        for(size_t i = 0; i < n->by_start.size(); i++) {
            printf("\t\t(%zu,%zu) ",
                   START(n->by_start[i]),
                   STOP(n->by_start[i]));
            printf("(%zu,%zu)\n",
                   START(n->by_stop[i]),
                   STOP(n->by_stop[i]));
        }
    }

    /* local print function */
    void print(IntervalTree<T>* n, int depth, char pos, int verbosity)
    {
        switch(verbosity) {
            case 1:
                print_verbosity_1(n, depth, pos);
                break;
            case 2:
                print_verbosity_2(n, depth, pos);
                break;
            case 3:
                print_verbosity_3(n, depth, pos);
                break;
            default:
                fprintf(stderr, "verbosity must be 1, 2, or 3\n");
                exit(EXIT_FAILURE);
        }
        depth++;
        if(LEFT(n) != nullptr) {
            print(LEFT(n), depth, 'l', verbosity);
        }
        if(RIGHT(n) != nullptr) {
            print(RIGHT(n), depth, 'r', verbosity);
        }
    }

protected:
    // the center position for this node
    long center = 0;
    // all intervals that overlap the center, sorted by start position
    std::vector<T*> by_start;
    // all intervals that overlap the center, sorted by stop position
    std::vector<T*> by_stop;
    // Child nodes
    IntervalTree* children[2] = {nullptr, nullptr};
    // Parent
    IntervalTree* parent = nullptr;
    // position relative to parent
    Orientation orientation = O_ROOT;

public:
    IntervalTree(
        std::vector<T*> intervals,
        IntervalTree<T>* new_parent = nullptr,
        Orientation new_orientation = O_ROOT
    ) :
        parent(new_parent), orientation(new_orientation)
    {

        center = get_center(intervals);

        std::vector<T*> left;
        std::vector<T*> right;

        /** iterate over intervals classifying and counting each */
        for(size_t i = 0; i < intervals.size(); i++) {
            T* v = intervals[i];
            switch(v->position_relative_to(center)) {
                case lo:
                    right.push_back(v);
                    break;
                case hi:
                    left.push_back(v);
                    break;
                case in:
                    by_stop.push_back(v);
                    by_start.push_back(v);
                    break;
                default:
                    break;
            }
        }

        if(!by_start.empty())
            std::sort(by_start.begin(), by_start.end(), cmp_start);

        if(!by_stop.empty())
            std::sort(by_stop.begin(), by_stop.end(), cmp_stop);

        if (left.empty()) {
            children[0] = nullptr;
        } else {
            children[0] = new IntervalTree<T>(left, this, O_LEFT);
        }

        if (right.empty()) {
            children[1] = nullptr;
        } else {
            children[1] = new IntervalTree<T>(right, this, O_RIGHT);
        }
    }

    ~IntervalTree()
    {
        delete children[0];
        delete children[1];
    }

    template <class U>
    long count_overlaps(U* inv)
    {
        return count_overlaps(inv, this, 0);
    }

    long count_overlaps(long pnt)
    {
        return count_overlaps(pnt, this, 0);
    }

    template <class U>
    IntervalResult<T>* get_overlaps(U* inv)
    {
        auto res = new IntervalResult<T>;
        res->tree = this;
        get_overlaps(inv, this, res);
        return res;
    }

    IntervalResult<T>* get_point_overlaps(long pnt)
    {
        IntervalResult<T>* res = new IntervalResult<T>;
        res->tree = this;
        get_overlaps(pnt, this, res);
        return res;
    }


    /* public wrapper for real print function */
    void print(int verbosity)
    {
        print(this, 0, 'c', verbosity);
    }

};

#endif
