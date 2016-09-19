#include "interval_set.h"

template<typename T>
T* IntervalSet<T>::front()
{
    return inv.front();
}

template<typename T>
T* IntervalSet<T>::back()
{
    return inv.back();
}

template<typename T>
void IntervalSet<T>::add(T* x)
{
    inv.push_back(x);
}

template<typename T>
void IntervalSet<T>::print()
{
    for(auto &x : inv) {
        x->print();
    }
}

template<typename T>
bool IntervalSet<T>::empty()
{
    return inv.empty();
}

template<typename T>
size_t IntervalSet<T>::size()
{
    return inv.size();
}

template<typename T>
void IntervalSet<T>::clear()
{
    inv.clear();
    delete tree;
}

template<typename T>
void IntervalSet<T>::sort()
{
    std::sort(inv, cmp_start);
}

template<typename T>
void IntervalSet<T>::sort(int sort_method)
{
    switch(sort_method){
        case 0:
            std::sort(inv, cmp_start);
            break;
        case 1:
            std::sort(inv, cmp_start_reverse);
            break;
        case 2:
            std::sort(inv, cmp_stop);
            break;
        case 3:
            std::sort(inv, cmp_stop_reverse);
            break;
        default:
            fprintf(stderr, "Sort method must be in {0,1,2,3}\n");
            exit(EXIT_FAILURE);
    }
}

template <class T>
bool IntervalSet<T>::cmp_start(T* a, T* b)
{
    return ( a->pos[0] < b->pos[0] );
}

template <class T>
bool IntervalSet<T>::cmp_stop(T* a, T* b)
{
    return ( a->pos[1] < b->pos[1] );
}

template <class T>
bool IntervalSet<T>::cmp_start_reverse(T* a, T* b)
{
    return ( a->pos[0] > b->pos[0] );
}

template <class T>
bool IntervalSet<T>::cmp_stop_reverse(T* a, T* b)
{
    return ( a->pos[1] > b->pos[1] );
}

template<typename T>
IntervalResult<T>* IntervalSet<T>::get_region(Bound& bound, bool get_flank_overlaps)
{
    tree = build_tree(tree);

    // Search itree
    IntervalResult<T>* res = tree->get_overlaps(&bound);

    // I want everything that overlaps the flanks if I am doing a Block search.
    // For ContiguousSet searches, I currently only want overlapping intervals.
    // TODO find a better solution
    if (get_flank_overlaps) {
        res = add_whatever_overlaps_flanks(res);
    }

    return res;
}

template<typename T>
T* IntervalSet<T>::add_whatever_overlaps_flanks(T* res)
{
    // itree returns the flanks for queries that overlap nothing. However, I
    // need all the intervals that overlap these flanks as well.
    T* tmp_a = NULL;
    T* tmp_b = NULL;

    if (res->inbetween) {
        // If inbetween, itree should have returned the two flanking blocks
        if (res->iv.size() == 2) {
            tmp_a = res->tree->get_overlaps(res->iv[0]);
            tmp_b = res->tree->get_overlaps(res->iv[1]);
            res->iv.clear();
            res->iv = tmp_a->iv;
            res->iv.insert(
                res->iv.end(),
                tmp_b->iv.begin(),
                tmp_b->iv.end()
            );
        } else {
            fprintf(stderr, "itree is broken, should return exactly 2 intervals for inbetween cases\n");
            exit(EXIT_FAILURE);
        }
    } else if (res->leftmost || res->rightmost) {
        if (res->iv.size() == 1) {
            tmp_a = res->tree->get_overlaps(res->iv[0]);
            res->iv = tmp_a->iv;
        } else {
            fprintf(stderr, "itree is broken, should return only 1 interval for left/rightmost cases\n");
            exit(EXIT_FAILURE);
        }
    }

    if (tmp_a != NULL)
        delete tmp_a;

    if (tmp_b != NULL)
        delete tmp_b;

    return res;
}

template <class T>
void IntervalSet<T>::build_tree()
{
    if (tree == NULL) {

        std::vector<T*> intervals;

        for (T* c = inv.front(); c != NULL; c = c->next()) {
            intervals.push_back(c);
        }

        tree = new IntervalTree<T>(intervals);
    }
}
