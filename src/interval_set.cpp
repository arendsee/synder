#include "interval_set.h"

template <class T>
bool Interval<T>::cmp_start(T* a, T* b)
{
    return ( a->pos[0] < b->pos[0] );
}

template <class T>
bool Interval<T>::cmp_stop(T* a, T* b)
{
    return ( a->pos[1] < b->pos[1] );
}
