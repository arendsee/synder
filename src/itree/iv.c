#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>

#include "iv.h"

IV * iv_init(size_t available){
    IV * iv = (IV *)malloc(sizeof(IV)); 
    iv->size = 0;
    iv->available = available;
    iv->v = (Interval *)malloc(available * sizeof(Interval));
    return(iv);
}

void iv_free(IV * self){
    if(self){
        if(self->v){
            free(self->v);
        }
        free(self);
    }
}

void iv_add (IV * self, Interval dat){
    if(self->size == self->available){
        self->v = (Interval *)realloc(self->v, self->size * 2 * sizeof(Interval));
        self->available *= 2;
    }
    self->v[self->size] = dat;
    self->size++;
}

void iv_join (IV * self, IA * other){
    if(self->size + other->size > self->available){
        self->available = 2 * (self->size + other->size);
        self->v = (Interval *)realloc(self->v, self->available * sizeof(Interval));
    }
    memcpy(self->v + self->size, other->v, other->size * sizeof(other->v[0]));
    self->size += other->size;
}
