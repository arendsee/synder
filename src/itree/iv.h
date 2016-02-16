#ifndef __IV_H__
#define __IV_H__

#define IV_INITIAL_SIZE 8

#include <stdlib.h>

#include "interval.h"
#include "ia.h"

/** An automatically resizing vector of Interval objects  */
typedef struct {
    size_t available; /* total allocated memory */
    size_t size;      /* memory used */
    Interval * v;
} IV;

IV * iv_init(size_t);

void iv_add (IV *, Interval);

void iv_join (IV *, IA *);

void iv_free(IV *);

#endif
