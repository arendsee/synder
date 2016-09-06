#include "ia.hpp"

IA * init_IA(){
    IA * ia = (IA *)malloc(sizeof(IA));
    ia->size = 0;
    ia->v = NULL;
    return(ia);
}

IA * init_set_IA(size_t size){
    IA * ia = (IA *)malloc(sizeof(IA));
    ia->size = size;
    if(size > 0){
        ia->v = (Interval *)malloc(size * sizeof(Interval));
    } else {
        ia->v = NULL;
    }
    return(ia);
}

void free_IA(IA * ia){
    if(ia->v != NULL)
        free(ia->v);
    if(ia != NULL)
        free(ia);
}
