#include "interval.h"

bool overlap(uint a1, uint a2, uint b1, uint b2){
    return a1 <= b2 && a2 >= b1;
}

uint anchor(uint x, Contig * contig){
    Block ** blks = contig->block;
    uint lo = 0;
    uint hi = contig->size - 1;
    uint i = hi / 2;
    while(true){
        if(i <= hi && x > blks[i]->stop){
            lo = i;
            i = (i + hi) / 2 + 1;
        }
        else if(i > lo && x < blks[i]->start){
            hi = i;
            if(i == lo + 1)
                break;
            i = lo + (i - lo) / 2;
        }
        else {
            break;
        }
    }
    return i;
}
