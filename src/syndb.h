#ifndef __SYNDB_H__
#define __SYNDB_H__

#include "ftypes.h"

Synmap * load_synmap(FILE *);

void free_synmap(Synmap *);

void print_block(Block *);

void print_contig(Contig *);

void print_genome(Genome *);

void print_synmap(Synmap *);

#endif
