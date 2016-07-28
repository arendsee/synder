#ifndef __IA_H__
#define __IA_H__

#include <stdlib.h>

#include "interval.h"

/** Interval array */
typedef struct {
    size_t size;
    Interval * v;
} IA;

/** initialize to {.size=0, .v=NULL} */
IA * init_IA();

/** initialize to a size, allocate memory to v */
IA * init_set_IA(size_t size);

void free_IA(IA * ia);

#endif
