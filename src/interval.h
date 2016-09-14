#ifndef __INTERVAL_H__
#define __INTERVAL_H__

#include <stdio.h>
#include <stdlib.h>

#define START(inv) inv->pos[0]
#define STOP(inv) inv->pos[1]


/** position of Interval or point A relative to B */
typedef enum {lo = 0, in = 1, hi = 2} Pos;

template <class T>
class Interval
{
public:
    long pos[2];

    /** find position of point A relative to interval B (see Pos) */
    Pos position_relative_to(long a);

    /** find position of this interval other */
    template <class U>
    Pos position_relative_to(U* b);

    /** Determine whether intervals */
    template <class U>
    bool overlap(U* b);

    /** Calculate the length of the overlap of two intervals */
    template <class U>
    long overlap_length(U* b);

};

#endif
