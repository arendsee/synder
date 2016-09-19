#include "linked_interval.h"

template<class T>
T* LinkedInterval<T>::prev()
{
    return cor[0];
}

template<class T>
T* LinkedInterval<T>::next()
{
    return cor[1];
}

template<class T>
T* LinkedInterval<T>::prev_adj()
{
    return adj[0];
}

template<class T>
T* LinkedInterval<T>::next_adj()
{
    return adj[1];
}

template<class T>
T* LinkedInterval<T>::corner(size_t i)
{
    return cor[i];
}

template<class T>
T*      LinkedInterval<T>::get_over()
{
    return over;
}

template<class T>
Contig* LinkedInterval<T>::get_parent()
{
    return parent;
}

template<class T>
char    LinkedInterval<T>::get_strand()
{
    return strand;
}

template<class T>
size_t  LinkedInterval<T>::get_grpid()
{
    return grpid;
}

template<class T>
size_t  LinkedInterval<T>::get_id()
{
    return id;
}

template<class T>
static void link_homologs(T* a, T* b){
    a->over = b;
    b->over = a;
}
