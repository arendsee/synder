#ifndef __LINKED_INTERVAL_H__
#define __LINKED_INTERVAL_H__

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

    virtual T* next();
    virtual T* prev();
    virtual T* next_bs();
    virtual T* prev_bs();
    virtual T* next_adj();
    virtual T* prev_adj();
}

#endif
