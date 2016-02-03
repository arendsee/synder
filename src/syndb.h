#ifndef __SYNDB_H__
#define __SYNDB_H__

#include <stdlib.h>
#include <stdio.h>

#define uint unsigned int

typedef struct {
    uint start;
    uint stop;
    uint oseqid;
    uint oblkid;
    uint linkid;
} Block;

typedef struct {
    char * name;
    size_t size;
    Block ** block;
} Contig;

typedef struct {
    char * name;
    size_t size;
    Contig ** contig;
} Genome;

typedef struct {
    Genome ** genome;
} Synmap;

Synmap * load_synmap(FILE *);

void free_synmap(Synmap *);

void print_block(Block *);
void print_contig(Contig *);
void print_genome(Genome *);
void print_synmap(Synmap *);

#endif
