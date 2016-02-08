#ifndef __SYNMAP_H__
#define __SYNMAP_H__

#include "global.h"
#include "genome.h"

typedef struct {
    Genome ** genome;
} Synmap;

Synmap * init_synmap();

void free_synmap(Synmap *);

Synmap * load_synmap(FILE *);

void print_synmap(Synmap *);

#endif
