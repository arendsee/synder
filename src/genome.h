#ifndef __GENOME_H__
#define __GENOME_H__

#include "global.h"
#include "contig.h"

typedef struct {
    char * name;
    size_t size;
    Contig ** contig;
} Genome;

Genome * init_genome(char *, size_t);

void free_genome(Genome *);

void print_genome(Genome *);

#endif
