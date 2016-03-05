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



// ------------------------------------------------------------------
// Test code
// ------------------------------------------------------------------

bool test_iv(){
    IV * iv = iv_init(2);
    Interval * a = init_interval(12,22);
    Interval * b = init_interval(23,33);
    Interval * c = init_interval(34,44);

    printf("iv: vector realloc when size is passed\n");
    iv_add(iv, *a);
    iv_add(iv, *b);
    assert(iv->available == 2);
    iv_add(iv, *c);
    assert(iv->available == 4);

    printf("iv: correct addition of elements\n");
    assert(iv->v[0].start == 12);
    assert(iv->v[1].start == 23);
    assert(iv->v[2].start == 34);

    free(a);
    free(b);
    free(c);
    iv_free(iv);

    return true;
}
