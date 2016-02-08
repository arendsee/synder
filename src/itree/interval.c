#include <stdlib.h>

#include "interval.h"

/* compare Intervals by stop */
int cmp_stop(const void *ap, const void *bp){
    Interval a = * (Interval *) ap;
    Interval b = * (Interval *) bp;
    return((a.stop > b.stop) - (b.stop > a.stop));
}

/* compare Intervals by start */
int cmp_start(const void *ap, const void *bp){
    Interval a = * (Interval *) ap;
    Interval b = * (Interval *) bp;
    return((a.start > b.start) - (b.start > a.start));
}

/* find position of point A relative to interval B (see Pos) */
Pos point_overlap(unsigned int a, Interval b){
    if(a < b.start){
        return lo;
    }
    else if(a > b.stop){
        return hi;
    }
    else{
        return in;
    }
}

/* find position of interval A relative to interval B (see Pos) */
Pos interval_overlap(Interval a, Interval b){
    if(a.stop < b.start){
        return lo;
    }
    else if(a.start > b.stop){
        return hi;
    }
    else{
        return in;
    }
}
