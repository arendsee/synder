#ifndef __CONTIG_H__
#define __CONTIG_H__

#include "global.h"
#include "block.h"
#include "result.h"

typedef struct {
    char * name;
    size_t size;
    struct IntervalTree * itree;
    Block ** block;
} Contig;

Contig * init_contig(char *, size_t);

void free_contig(Contig *);

void print_contig(Contig *);

/**
 * Find index of downstream Block nearest the query point
 */
uint anchor(uint, Contig *);

/**
 * Given two points, find all blocks overlapping them
 */
Contig * get_overlapping(uint, uint, Contig *);

/**
 * Given two points, find the blocks flanking them
 */
Contig * get_flanks(uint, uint, Contig *, uint, uint);

/**
 * Given two points, find the number of blocks they overlap
 */
uint count_overlaps(uint, uint, Contig *);

/**
 * Given two points, get the expected location on the target
 */
Result map(uint, uint, Contig *);

#endif
