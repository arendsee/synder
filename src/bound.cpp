#include "bound.h"

Bound::Bound(){
    pos[0] = 0;
    pos[1] = 0;
}

Bound::Bound(long start, long stop){
    pos[0] = start;
    pos[1] = stop;
}

Bound::~Bound(){ }
