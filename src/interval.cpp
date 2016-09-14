#include "interval.h"

template <class T>
Pos Interval<T>::position_relative_to(long a)
{
    if (a < pos[0]) {
        return lo;
    } else if (a > pos[1]) {
        return hi;
    } else {
        return in;
    }
}

template <class T>
template <class U>
Pos Interval<T>::position_relative_to(U* b)
{
    if (pos[1] < b->pos[0]) {
        return lo;
    } else if (pos[0] > b->pos[1]) {
        return hi;
    } else {
        return in;
    }
}

template <class T>
template <class U>
long Interval<T>::overlap_length(U* b)
{
    long a1 = pos[0];
    long a2 = pos[1];
    long b1 = b->pos[0];
    long b2 = b->pos[1];

    // If the intervals overlap
    if(a1 <= b2 && b1 <= a2) {
        // Find the lower bound of the overlapping region
        long a = a1 > b1 ? a1 : b1;
        // Find the upper bound of the overlapping region
        long b = a2 > b2 ? b2 : a2;
        // Return the overlapping interval length
        return b - a + 1;
    } else {
        return 0;
    }
}

template <class T>
template <class U>
bool Interval<T>::overlap(U* b)
{
    long a1 = pos[0];
    long a2 = pos[1];
    long b1 = b->pos[0];
    long b2 = b->pos[1];

    return (a1 <= b2) && (a2 >= b1);
}
