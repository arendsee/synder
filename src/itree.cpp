#include "itree.h"

IntervalTree::IntervalTree(
    Interval * new_interval_pool,
    size_t new_size
)
    : size(new_size), interval_pool(new_interval_pool)
{
    std::vector<Interval*> intervals(size);
    for (size_t i = 0; i < size; i++)
    {
        intervals[i] = &interval_pool[i];
    }
    root = new IntervalTreeNode(intervals);
}

IntervalTree::~IntervalTree()
{
    if (interval_pool != NULL)
        free(interval_pool);
    if (root != NULL)
        delete root;
}
 
void IntervalTree::print(int verbosity)
{
    root->print(verbosity);
}

long IntervalTree::count_overlaps(long pnt)
{
    return count_overlaps(pnt, root, 0);
}

long IntervalTree::count_overlaps(Interval * inv)
{
    return count_overlaps(inv, root, 0);
}

IntervalResult * IntervalTree::get_overlaps(long pnt)
{
    IntervalResult * res = new IntervalResult;
    get_overlaps(pnt, root, res);
    return res;
}

IntervalResult * IntervalTree::get_overlaps(Interval * inv)
{
    IntervalResult * res = new IntervalResult;
    get_overlaps(inv, root, res);
    return res;
}


long IntervalTree::count_overlaps(long pnt, IntervalTreeNode * tree, long count)
{
    if (pnt >= tree->center)
    {
        for (long long i = T_SIZE(tree) - 1; i >= 0 ; i--)
        {
            if (pnt <= T_STOP_STOP(tree, i))
            {
                count++;
            }
            else
            {
                break;
            }
        }
        if (RIGHT(tree) != NULL)
            return count_overlaps(pnt, RIGHT(tree), count);
    }
    else
    {
        for (size_t i = 0; i < T_SIZE(tree); i++)
        {
            if (pnt >= T_START_START(tree, i))
            {
                count++;
            }
            else
            {
                break;
            }
        }
        if (LEFT(tree) != NULL)
            return count_overlaps(pnt, LEFT(tree), count);
    }
    return count;
}

long IntervalTree::count_overlaps(Interval * inv, IntervalTreeNode * tree, long count)
{
    if (tree == NULL)
        return count;
    switch (inv->overlap(tree->center))
    {
    case lo:
        for (long long i = T_SIZE(tree) - 1; i >= 0 ; i--)
        {
            if (inv->start <= T_STOP_STOP(tree, i))
            {
                count++;
            }
            else
            {
                break;
            }
        }
        return count_overlaps(inv, RIGHT(tree), count);
    case hi:
        for (size_t i = 0; i < T_SIZE(tree); i++)
        {
            if (inv->stop >= T_START_START(tree, i))
            {
                count++;
            }
            else
            {
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

void IntervalTree::set_nearest_opposing_interval(IntervalTreeNode * node, IntervalResult * results, Pos pos)
{

    Orientation orientation = node->orientation;

    if ((pos == hi && orientation == O_LEFT) ||
            (pos == lo && orientation == O_RIGHT))
    {
        node = node->parent;
    }
    else
    {
        while (node->orientation != O_ROOT && node->orientation == orientation)
        {
            node = node->parent;
        }
        node = node->parent;
    }

    if (node == NULL)
    {
        if (pos == lo)
        {
            results->leftmost = true;
        }
        else
        {
            results->rightmost = true;
        }
    }
    else
    {
        results->inbetween = true;
        if (orientation == O_RIGHT)
        {
            results->iv.push_back(node->by_start.front());
        }
        else if (orientation == O_LEFT)
        {
            results->iv.push_back(node->by_stop[node->by_stop.size() - 1]);
        }
    }
}

void IntervalTree::get_overlaps(long pnt, IntervalTreeNode * tree, IntervalResult * results)
{
    if (pnt >= tree->center)
    {
        for (long long i = T_SIZE(tree) - 1; i >= 0 ; i--)
        {
            if (pnt <= T_STOP_STOP(tree, i))
            {
                results->iv.push_back(T_STOP(tree, i));
            }
            else
            {
                break;
            }
        }
        if (RIGHT(tree) != NULL)
        {
            get_overlaps(pnt, RIGHT(tree), results);
        }
        else if (R_SIZE(results) != 0)
        {
            results->iv.push_back(LAST_STOP(tree));
            set_nearest_opposing_interval(tree, results, hi);
        }
    }
    else
    {
        for (size_t i = 0; i < T_SIZE(tree); i++)
        {
            if (pnt >= T_START_START(tree, i))
            {
                results->iv.push_back(T_START(tree, i));
            }
            else
            {
                break;
            }
        }
        if (LEFT(tree) != NULL)
        {
            get_overlaps(pnt, LEFT(tree), results);
        }
        else if (R_SIZE(results) != 0)
        {
            results->iv.push_back(FIRST_START(tree));
            set_nearest_opposing_interval(tree, results, lo);
        }
    }
}

IntervalResult * IntervalTree::get_overlaps(Interval * inv, IntervalTreeNode * tree, IntervalResult * results)
{
    if (tree == NULL)
        return results;
    switch (inv->overlap(tree->center))
    {
    case lo: // center lower than interval start
        for (long long i = T_SIZE(tree) - 1; i >= 0 ; i--)
        {
            if (inv->start <= T_STOP_STOP(tree, i))
            {
                results->iv.push_back(T_STOP(tree, i));
            }
            else
            {
                break;
            }
        }
        // Reach a leaf and still have no overlaps
        if (RIGHT(tree) == NULL && results->iv.size() == 0)
        {
            results->iv.push_back(LAST_STOP(tree));
            // get nearest interval on the other side
            set_nearest_opposing_interval(tree, results, hi);
        }
        return get_overlaps(inv, RIGHT(tree), results);
    case hi:
        for (size_t i = 0; i < T_SIZE(tree); i++)
        {
            if (inv->stop >= T_START_START(tree, i))
            {
                results->iv.push_back(tree->by_start[i]);
            }
            else
            {
                break;
            }
        }
        // Reach a leaf and still have no overlaps
        if (LEFT(tree) == NULL && results->iv.size() == 0)
        {
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
