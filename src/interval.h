#ifndef __INTERVAL_H__
#define __INTERVAL_H__

#include <stdio.h>
#include <stdlib.h>

/** position of Interval or point A relative to B */
typedef enum {lo = 0, in = 1, hi = 2} Pos;

class Interval
{
public:
    long start;
    long stop;
    void * link; /* a pointer to arbitrary related data */

    Interval(long start, long stop);

    void print();

    /** find position of point A relative to interval B (see Pos) */
    Pos overlap(long A);

    /** find position of this interval other */
    Pos overlap(Interval * other);

    /** compare intervals by stop */
    static bool cmp_stop(Interval * a, Interval * b);

    /** compare intervals by start */
    static bool cmp_start(Interval * a, Interval * b);

};

#endif
