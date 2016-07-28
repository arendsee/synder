#ifndef __IV_H__
#define __IV_H__

#define IV_INITIAL_SIZE 8

#include <stdlib.h>
#include <stdbool.h>


#include "interval.h"
#include "ia.h"

/** An automatically resizing vector of Interval objects  */
typedef struct {
    size_t available; /* total allocated memory */
    size_t size;      /* memory used */
    Interval * v;
} IV;

IV * init_IV(size_t);

void print_IV(IV * self);

void add_IV(IV *, Interval);

void join_IV(IV *, IA *);

void free_IV(IV *);

#endif
