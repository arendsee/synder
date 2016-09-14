#ifndef __INTERVAL_SET_H__
#define __INTERVAL_SET_H__

template <class T>
class IntervalSet
{
public:
    std::vector<T*> inv;
    IntervalTree<T>* itree;

    virtual T* begin();
    virtual void add(T* );
    virtual void print();

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
    IntervalResult<T>* get_region(long a, long b);

    /** Sort blocks */
    void sort(bool by_stop);
    void sort();

    /** compare intervals by stop - used in std::sort */
    static bool cmp_stop(T* a, T* b);

    /** compare intervals by start - used in std::sort */
    static bool cmp_start(T* a, T* b);

}

#endif
