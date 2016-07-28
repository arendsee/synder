#include <stdlib.h>
#include <stdio.h>

#include "interval.h"

Interval * init_Interval(uint start, uint stop){
    Interval * inv = (Interval*)malloc(sizeof(Interval));
    inv->start = start;
    inv->stop = stop;
    inv->link = NULL;
    return(inv);
}

void print_Interval(Interval * interval){
    printf("%u %u\n", interval->start, interval->stop);
}

int cmp_stop(const void *ap, const void *bp){
    Interval a = * (Interval *) ap;
    Interval b = * (Interval *) bp;
    return((a.stop > b.stop) - (b.stop > a.stop));
}

int cmp_start(const void *ap, const void *bp){
    Interval a = * (Interval *) ap;
    Interval b = * (Interval *) bp;
    return((a.start > b.start) - (b.start > a.start));
}

Pos point_overlap(uint a, Interval * b){
    if(a < b->start){
        return lo;
    }
    else if(a > b->stop){
        return hi;
    }
    else{
        return in;
    }
}

Pos interval_overlap(Interval * a, Interval * b){
    if(a->stop < b->start){
        return lo;
    }
    else if(a->start > b->stop){
        return hi;
    }
    else{
        return in;
    }
}
