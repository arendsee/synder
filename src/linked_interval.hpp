#ifndef __LINKED_INTERVAL_H__
#define __LINKED_INTERVAL_H__

#include "global.h"
#include "feature.h"
#include "bound.h"

#include <array>


template<class T> class LinkedInterval
{
protected:
    // next and prev elements by start and stop
    std::array<T*, 4> cor = {{ nullptr }};

    // adjacent non-overlapping block
    std::array<T*, 2> adj = {{ nullptr }};

public:
    // homologous element
    T* over = nullptr;

    // parent of this interval
    Feature* parent = nullptr;

    // score provided by synteny program
    double score = 0;

    // strand relative to query
    char strand = '+';

    // overlapping group id
    long grpid = 0;

    // unique id for element
    size_t id = 0;

    LinkedInterval() { }

    LinkedInterval(
        Feature* t_parent,
        double   t_score,
        char     t_strand,
        size_t   t_linkid
    )
        :
        parent ( t_parent ),
        score  ( t_score  ),
        strand ( t_strand ),
        id     ( t_linkid )
    { }

    ~LinkedInterval() { }

    T* prev()        { return cor[0]; }
    T* next()        { return cor[1]; }
    T* prev_bystop() { return cor[2]; }
    T* next_bystop() { return cor[3]; }
    T* prev_adj()    { return adj[0]; }
    T* next_adj()    { return adj[1]; }

    T* corner(size_t i)
    {
        try {
            return cor.at(i);
        } catch (const std::out_of_range& e) {
            std::cerr << "Index error in " << __func__ << std::endl;
            return nullptr;
        }
    }

    T* corner_adj(size_t i)
    {
        try {
            return adj.at(i);
        } catch (const std::out_of_range& e) {
            std::cerr << "Index error in " << __func__ << std::endl;
            return nullptr;
        }
    }

    static void link_homologs(T* a, T* b)
    {
        a->over = b;
        b->over = a;
    }
};

#endif
