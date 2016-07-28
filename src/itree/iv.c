#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>

#include "iv.h"

IV * init_IV(size_t available){
    IV * iv = (IV *)malloc(sizeof(IV)); 
    if(iv == NULL){
        exit(EXIT_FAILURE);
    }
    iv->size = 0;
    iv->available = available;
    iv->v = (Interval *)malloc(available * sizeof(Interval));
    if(iv->v == NULL){
        exit(EXIT_FAILURE);
    }
    return(iv);
}

void free_IV(IV * self){
    if(self != NULL){
        if(self->v != NULL){
            free(self->v);
        }
        free(self);
    }
}

void print_IV(IV * self){
    for(size_t i = 0; i < self->size; i++){
        printf("(%u, %u) ", self->v[i].start, self->v[i].stop); 
    }
    printf("\n");
}

void add_IV (IV * self, Interval dat){
    if(self->size == self->available){
        self->v = (Interval *)realloc(self->v, self->size * 2 * sizeof(Interval));
        if(self->v == NULL){
            exit(EXIT_FAILURE);
        }
        self->available *= 2;
    }
    self->v[self->size] = dat;
    self->size++;
}

void join_IV (IV * self, IA * other){
    if(self->size + other->size > self->available){
        self->available = 2 * (self->size + other->size);
        self->v = (Interval *)realloc(self->v, self->available * sizeof(Interval));
        if(self->v == NULL){
            exit(EXIT_FAILURE);
        }
    }
    memcpy(self->v + self->size, other->v, other->size * sizeof(other->v[0]));
    self->size += other->size;
}
