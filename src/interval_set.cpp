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


IntervalResult<T>* get_region(long a, long b){
    // TODO implement
}

T* begin(){
    return new T;
    // TODO implement
}

void add(T* ){
    // TODO implement
}

void print(){
    // TODO implement
}

void sort(bool by_stop){
    // TODO implement
}

void sort(){
    // TODO implement
}
