#ifndef __LINKED_INTERVAL_H__
#define __LINKED_INTERVAL_H__

#include "global.h"
#include "feature.h"

#include <array>


template<class T> class LinkedInterval
{
protected:
    // homologous element
    T* over = nullptr;

    // parent of this interval
    Feature* parent = nullptr;

    // strand relative to query
    char strand = '+';

    // next and prev elements by start and stop
    std::array<T*, 4> cor = {nullptr, nullptr, nullptr, nullptr};

    // adjacent non-overlapping block
    std::array<T*, 2> adj;

    // score provided by synteny program
    double score = 0;

    // overlapping group id
    size_t grpid = 0;

    // unique id for element
    size_t id = 0;

public:
    LinkedInterval() {}

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
    {}

    T* prev()
    {
        return cor[0];
    }

    T* next()
    {
        return cor[1];
    }

    T* prev_adj()
    {
        return adj[0];
    }

    T* next_adj()
    {
        return adj[1];
    }

    T* corner(size_t i)
    {
        try {
            return cor.at(i);
        } catch (const std::out_of_range& e) {
            cerr << "Index error in " << __func__ << endl;
            return NULL;
        }
    }

    T* get_over()
    {
        return over;
    }

    Feature* get_parent()
    {
        return parent;
    }

    char get_strand()
    {
        return strand;
    }

    size_t get_grpid()
    {
        return grpid;
    }

    size_t get_id()
    {
        return id;
    }

    static void link_homologs(T* a, T* b)
    {
        a->over = b;
        b->over = a;
    }
};

#endif
