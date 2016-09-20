#ifndef __INTERVAL_H__
#define __INTERVAL_H__

#include <stdio.h>
#include <stdlib.h>
#include <array>

#define START(inv) inv->pos[0]
#define STOP(inv) inv->pos[1]

/** position of Interval or point A relative to B */
typedef enum {lo = 0, in = 1, hi = 2} Pos;

template <class T>
class Interval
{
public:

    std::array<long, 2> pos = { 0, 0 };

    Interval() {}

    Interval(long t_start, long t_stop)
        :
        pos ( { t_start, t_stop }  )
    {}

    /** find position of point A relative to interval B (see Pos) */
    Pos position_relative_to(long a)
    {
        if (a < pos[0]) {
            return lo;
        } else if (a > pos[1]) {
            return hi;
        } else {
            return in;
        }
    }

    /** find position of this interval other */
    template <class U>
    Pos position_relative_to(U* other)
    {
        if (pos[1] < other->pos[0]) {
            return lo;
        } else if (pos[0] > other->pos[1]) {
            return hi;
        } else {
            return in;
        }
    }

    /** Determine whether intervals */
    template <class U>
    bool overlap(U* other)
    {
        long a1 = pos[0];
        long a2 = pos[1];
        long b1 = other->pos[0];
        long b2 = other->pos[1];

        // If the intervals overlap
        if(a1 <= b2 && b1 <= a2) {
            // Find the lower bound of the overlapping region
            long a = a1 > b1 ? a1 : b1;
            // Find the upper bound of the overlapping region
            long b = a2 > b2 ? b2 : a2;
            // Return the overlapping interval length
            return b - a + 1;
        } else {
            return 0;
        }
    }

    /** Calculate the length of the overlap of two intervals */
    template <class U>
    long overlap_length(U* other)
    {
        long a1 = pos[0];
        long a2 = pos[1];
        long b1 = other->pos[0];
        long b2 = other->pos[1];

        return (a1 <= b2) && (a2 >= b1);
    }

};

#endif
