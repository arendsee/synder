#ifndef __LINKED_INTERVAL_H__
#define __LINKED_INTERVAL_H__

#include "global.h"

template<class T>
class LinkedInterval
{
private:
    T*      over;   // homologous element
    Contig* parent; // Contig parent of this interval
    char    strand; // strand relative to query
    T*      cor[4]; // next and prev elements by start and stop
    T*      adj[2]; // adjacent non-overlapping block
    double  score;  // score provided by synteny program
    size_t  grpid;  // overlapping group id;
    size_t  id;     // unique id for element

public:

    T* prev()         { return cor[0]; }
    T* next()         { return cor[1]; }
    T* prev_adj()     { return adj[1]; }
    T* next_adj()     { return adj[1]; }
    T* cor(size_t i)  { return cor[i]; }

    T* over()         { return over;   }
    Conitig* parent() { return parent; }
    char strand       { return strand; }
    size_t grpid      { return grpid;  }
    size_t id         { return id;     }

};

#endif
