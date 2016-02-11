#ifndef __IV_H__
#define __IV_H__

#include <stdlib.h>

#include "interval.h"

/** An automatically resizing vector of Interval objects  */
typedef struct {
    size_t available; /* total allocated memory */
    size_t size;      /* memory used */
    Interval * data;
} IV;

IV * iv_init(size_t);

void iv_add (IV *, Interval);

void iv_free(IV *);

#endif
