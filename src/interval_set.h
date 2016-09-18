#ifndef __INTERVAL_SET_H__
#define __INTERVAL_SET_H__

template <class T>
class IntervalSet
{
private:
    std::vector<T*> inv;
    IntervalTree<T>* tree;

    /** interval comparators - used in std::sort */
    bool cmp_start(T* a, T* b);
    bool cmp_stop(T* a, T* b);
    bool cmp_start_reverse(T* a, T* b);
    bool cmp_stop_reverse(T* a, T* b);

public:

    // All of these are simple wrappers for std:vector inv
    virtual T*     front();
    virtual T*     back();
    virtual bool   empty();
    virtual size_t size();
    virtual void   clear();

    // Iterate through intervals, calling print on each
    virtual void print();

    // wrapper for std::vector.push_back(T*)
    virtual void add(T*);

    /** Sort inv by pos[0]*/
    void sort();

    /** Sort inv in various ways
     *
     * char sort_method can be
     * - 0 = sort forward by start (same as sort())
     * - 1 = sort reverse by start
     * - 2 = sort forward by stop
     * - 3 = sort reverse by stop
     *
     */
    void sort(int sort_method);

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
    IntervalResult<T>* get_region(Bound& b);
};

#endif
