#include <stdlib.h>

#include "iv.h"

IV * iv_init(size_t available){
    IV * v = (IV *)malloc(sizeof(IV)); 
    v->size = 0;
    v->available = available;
    v->data = (Interval *)malloc(available * sizeof(Interval));
    return(v);
}

void iv_free(IV * self){
    if(self){
        if(self->data){
            free(self->data);
        }
        free(self);
    }
}

void iv_add (IV * self, Interval dat){
    if(self->size == self->available){
        self->data = (Interval *)realloc(self->data, self->size * 2 * sizeof(Interval));
        self->available *= 2;
    }
    self->data[self->size] = dat;
    self->size++;
}
