#include <stdlib.h>

#include "ia.h"
#include "interval.h"

IA * init_ia(){
    IA * ia = (IA *)malloc(sizeof(IA));
    ia->size = 0;
    ia->v = NULL;
    return(ia);
}

IA * init_set_ia(size_t size){
    IA * ia = (IA *)malloc(sizeof(IA));
    ia->size = size;
    if(size > 0){
        ia->v = (Interval *)malloc(size * sizeof(Interval));
    } else {
        ia->v = NULL;
    }
    return(ia);
}

void free_ia(IA * ia){
    if(ia->v)
        free(ia->v);
    if(ia)
        free(ia);
}
