#ifndef __LINKED_INTERVAL_H__
#define __LINKED_INTERVAL_H__

#include "global.h"

// Forward declarations
class Contig;

template<class T> class LinkedInterval
{
// TODO remove friendship
friend class ManyBlocks;

protected:
    T*      over;   // homologous element
    Contig* parent; // Contig parent of this interval
    char    strand; // strand relative to query
    T*      cor[4]; // next and prev elements by start and stop
    T*      adj[2]; // adjacent non-overlapping block
    double  score;  // score provided by synteny program
    size_t  grpid;  // overlapping group id;
    size_t  id;     // unique id for element

public:
    T* prev();
    T* next();
    T* prev_adj();
    T* next_adj();
    T* corner(size_t i);

    T*      get_over();
    Contig* get_parent();
    char    get_strand();
    size_t  get_grpid();
    size_t  get_id();

    static void link_homologs(T* a, T* b);
};

#endif
