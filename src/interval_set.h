#ifndef __INTERVAL_SET_H__
#define __INTERVAL_SET_H__

#include <vector>

#include "global.h"
#include "interval_tree.h"

/** A container for LinkedIntervals */
template <class T>
class IntervalSet
{
private:
    IntervalResult<T>* add_whatever_overlaps_flanks(IntervalResult<T>* res)
    {
        // itree returns the flanks for queries that overlap nothing. However, I
        // need all the intervals that overlap these flanks as well.
        IntervalResult<T>* tmp_a = nullptr;
        IntervalResult<T>* tmp_b = nullptr;

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

        delete tmp_a;
        delete tmp_b;

        return res;
    }

    void build_tree()
    {
        if (tree == nullptr) {
            tree = new IntervalTree<T>(inv);
        }
    }

protected:
    std::vector<T*> inv;
    IntervalTree<T>* tree;
    Feature* parent;

    static bool cmp_start         (T* a, T* b) { return ( a->pos[0] < b->pos[0] ); }
    static bool cmp_stop          (T* a, T* b) { return ( a->pos[1] < b->pos[1] ); }
    static bool cmp_start_reverse (T* a, T* b) { return ( a->pos[0] > b->pos[0] ); }
    static bool cmp_stop_reverse  (T* a, T* b) { return ( a->pos[1] > b->pos[1] ); }

public:
    virtual ~IntervalSet()
    {
        delete tree;
    }

    virtual T*     front() { return inv.front(); }
    virtual T*     back()  { return inv.back();  }
    virtual bool   empty() { return inv.empty(); }
    virtual size_t size()  { return inv.size();  }

    virtual void clear()
    {
        inv.clear();
        delete tree;
    }

    // Iterate through intervals, calling print on each
    virtual void print()
    {
        for(auto &x : inv) {
            x->print();
        }
    }

    // wrapper for std::vector.push_back(T*)
    virtual void add(T* x)
    {
        inv.push_back(x);
    }

    /** Sort inv by pos[0]*/
    void sort()
    {
        std::sort(inv.begin(), inv.end(), IntervalSet<T>::cmp_start);
    }

    /** Sort inv in various ways
     *
     * char sort_method can be
     * - 0 = sort forward by start (same as sort())
     * - 1 = sort reverse by start
     * - 2 = sort forward by stop
     * - 3 = sort reverse by stop
     *
     */
    void sort(int sort_method)
    {
        switch(sort_method) {
            case 0:
                std::sort(inv.begin(), inv.end(), IntervalSet<T>::cmp_start);
                break;
            case 1:
                std::sort(inv.begin(), inv.end(), IntervalSet<T>::cmp_start_reverse);
                break;
            case 2:
                std::sort(inv.begin(), inv.end(), IntervalSet<T>::cmp_stop);
                break;
            case 3:
                std::sort(inv.begin(), inv.end(), IntervalSet<T>::cmp_stop_reverse);
                break;
            default:
                throw "Sort method must be in {0,1,2,3}";
        }
    }

    /** Given two points, find all blocks overlapping or flanking them
     *
     * If there is at least one overlapping block, return all overlapping blocks
     *
     * Otherwise, return the blocks above and below the input region
     *
     * If there is only one flanking interval (i.e., the query is beyond any
     * syntenic interval), return just the nearest interval.
     *
     */
    template<class U>
    IntervalResult<T>* get_region(U& bound, bool get_flank_overlaps)
    {
        build_tree();

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

    template<class U>
    long count_overlaps(U* u)
    {
        build_tree();
        return tree->count_overlaps(u);
    }

    long count_point_overlaps(long pnt)
    {
        build_tree();
        return tree->count_point_overlaps(pnt);
    }
};

#endif
