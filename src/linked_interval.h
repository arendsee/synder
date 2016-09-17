#ifndef __LINKED_INTERVAL_H__
#define __LINKED_INTERVAL_H__

#include "global.h"

template<class T>
class LinkedInterval
{
public:
    T*      over;   // homologous element
    Contig* parent; // Contig parent of this interval
    char    strand; // strand relative to query
    T*      cor[4]; // next and prev elements by start and stop
    T*      adj[2]; // adjacent non-overlapping block
    double  score;  // score provided by synteny program
    size_t  grpid;  // overlapping group id;
    size_t  id;     // unique id for element
};

#endif
