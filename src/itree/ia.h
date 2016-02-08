#ifndef __IA_H__
#define __IA_H__

#include <stdlib.h>

#include "interval.h"

/** Interval array */
typedef struct {
    size_t size;
    Interval * v;
} IA;

/* initialize to {.size=0, .v=NULL} */
IA * init_ia();

/* initialize to a size, allocate memory to v */
IA * init_set_ia(size_t size);

void free_ia(IA * ia);

#endif
