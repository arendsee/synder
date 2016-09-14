#include "itree.h"

template <class T>
IntervalTree<T>::IntervalTree(
    std::vector<T*> intervals,
    IntervalTree<T>* new_parent,
    Orientation new_orientation
) :
    parent(new_parent), orientation(new_orientation)
{

    center = get_center(intervals);

    std::vector<T*> left;
    std::vector<T*> right;

    /* iterate over intervals classifying and counting each */
    T* v;
    for(size_t i = 0; i < intervals.size(); i++) {
        v = intervals[i];
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
        std::sort(by_start.begin(), by_start.end(), Interval<T>::cmp_start);

    if(!by_stop.empty())
        std::sort(by_stop.begin(), by_stop.end(), Interval<T>::cmp_stop);

    if (left.empty()) {
        children[0] = NULL;
    } else {
        children[0] = new IntervalTree<T>(left, this, O_LEFT);
    }

    if (right.empty()) {
        children[1] = NULL;
    } else {
        children[1] = new IntervalTree<T>(right, this, O_RIGHT);
    }
}


template <class T>
IntervalTree<T>::~IntervalTree()
{
    if (children[0] != NULL)
        delete children[0];

    if (children[1] != NULL)
        delete children[1];
}


template <class T>
long IntervalTree<T>::count_overlaps(long pnt)
{
    return count_overlaps(pnt, this, 0);
}


template <class T>
long IntervalTree<T>::count_overlaps(Bound* inv)
{
    return count_overlaps(inv, this, 0);
}


template <class T>
IntervalResult<T>* IntervalTree<T>::get_overlaps(long pnt)
{
    IntervalResult<T>* res = new IntervalResult<T>;
    res->tree = this;
    get_overlaps(pnt, this, res);
    return res;
}


template <class T>
IntervalResult<T>* IntervalTree<T>::get_overlaps(Bound* inv)
{
    auto res = new IntervalResult<T>;
    res->tree = this;
    get_overlaps(inv, this, res);
    return res;
}


template <class T>
long IntervalTree<T>::count_overlaps(long pnt, IntervalTree<T>* tree, long count)
{
    if (pnt >= tree->center) {
        for (long long i = T_SIZE(tree) - 1; i >= 0 ; i--) {
            if (pnt <= T_STOP_STOP(tree, i)) {
                count++;
            } else {
                break;
            }
        }
        if (RIGHT(tree) != NULL)
            return count_overlaps(pnt, RIGHT(tree), count);
    } else {
        for (size_t i = 0; i < T_SIZE(tree); i++) {
            if (pnt >= T_START_START(tree, i)) {
                count++;
            } else {
                break;
            }
        }
        if (LEFT(tree) != NULL)
            return count_overlaps(pnt, LEFT(tree), count);
    }
    return count;
}


template <class T>
long IntervalTree<T>::count_overlaps(Bound* inv, IntervalTree<T>* tree, long count)
{
    if (tree == NULL)
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


template <class T>
void IntervalTree<T>::set_nearest_opposing_interval(IntervalTree<T>* node, IntervalResult<T>* results, Pos pos)
{

    Orientation orientation = node->orientation;

    if ((pos == hi && orientation == O_LEFT) ||
            (pos == lo && orientation == O_RIGHT)) {
        node = node->parent;
    } else {
        while (node->orientation != O_ROOT && node->orientation == orientation) {
            node = node->parent;
        }
        node = node->parent;
    }

    if (node == NULL) {
        if (pos == lo) {
            results->leftmost = true;
        } else {
            results->rightmost = true;
        }
    } else {
        results->inbetween = true;
        if (orientation == O_RIGHT) {
            results->iv.push_back(node->by_start.front());
        } else if (orientation == O_LEFT) {
            results->iv.push_back(node->by_stop[node->by_stop.size() - 1]);
        }
    }
}


template <class T>
void IntervalTree<T>::get_overlaps(long pnt, IntervalTree<T>* tree, IntervalResult<T>* results)
{
    if (pnt >= tree->center) {
        for (long long i = T_SIZE(tree) - 1; i >= 0 ; i--) {
            if (pnt <= T_STOP_STOP(tree, i)) {
                results->iv.push_back(T_STOP(tree, i));
            } else {
                break;
            }
        }
        if (RIGHT(tree) != NULL) {
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
        if (LEFT(tree) != NULL) {
            get_overlaps(pnt, LEFT(tree), results);
        } else if (R_SIZE(results) != 0) {
            results->iv.push_back(FIRST_START(tree));
            set_nearest_opposing_interval(tree, results, lo);
        }
    }
}

template <class T>
IntervalResult<T>* IntervalTree<T>::get_overlaps(Bound* inv, IntervalTree<T>* tree, IntervalResult<T>* results)
{
    if (tree == NULL)
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
            if (RIGHT(tree) == NULL && results->iv.size() == 0) {
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
            if (LEFT(tree) == NULL && results->iv.size() == 0) {
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


/*
 * Select a point at the center of the middle interval.
 * This guarantees at least one interval overlaps each node.
 * If the intervals are sorted, it also favors (but doesn't guarantee) a
 * balanced tree.
 */
template <class T>
long IntervalTree<T>::get_center(std::vector<Bound*> intr)
{
    // get the central index
    long i = intr.size() / 2;
    // get the center point on this index
    long x = (STOP(intr[i]) - START(intr[i])) / 2 + START(intr[i]);
    return x;
}


/* write tree and center */
template <class T>
void IntervalTree<T>::print_verbosity_1(IntervalTree<T>* n, int depth, char pos)
{
    printf("%*d - %c%zu\n", depth * 2, depth, pos, n->center);
}


/* write tree, center, and start-sorted */
template <class T>
void IntervalTree<T>::print_verbosity_2(IntervalTree<T>* n, int depth, char pos)
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
template <class T>
void IntervalTree<T>::print_verbosity_3(IntervalTree<T>* n, int depth, char pos)
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
template <class T>
void IntervalTree<T>::print(IntervalTree<T>* n, int depth, char pos, int verbosity)
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
    if(LEFT(n) != NULL) {
        print(LEFT(n), depth, 'l', verbosity);
    }
    if(RIGHT(n) != NULL) {
        print(RIGHT(n), depth, 'r', verbosity);
    }
}


/* public wrapper for real print function */
template <class T>
void IntervalTree<T>::print(int verbosity)
{
    print(this, 0, 'c', verbosity);
}
