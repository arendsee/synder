#ifndef __CONTIG_H__
#define __CONTIG_H__

#include <stdlib.h>
#include <stdbool.h>

#include "block.h"

#ifndef uint
#define uint unsigned int
#endif

/** Contiguous sequence object containing list of Block structures and an
 * interval tree to search them.
 *
 * Fields:
 * - name  - a unique name for this contig (e.g. "Chr1")
 * - size  - number of blocks in the block array
 * - block - array of pointers to Block objects
 * - itree - an IntervalTree that allows log(n)+m search of overlapping intervals
 *
 * */
typedef struct {
    char * name;
    struct IntervalTree * itree;
    size_t size;
    Block ** block;
    Block ** by_stop;
    uint start_sorted:1;
    uint stop_sorted:1;
} Contig;

/** Allocate memory for a contig and set each field.
 *
 * The itree field, which can hold an interval tree for (log(n)+m) overlap
 * searching, is initialized to NULL. Building this tree is expensive (n
 * log(n)) so will not be done unless needed.
 *
 * Memory is allocated for pointers to *size* blocks, but they are not
 * initialized.
 *
 * @param name name of this Contig (e.g. "Chr1")
 * @param size number of Block structs this Contig will hold
 *
 * @return pointer to a new Contig
 *
 * */
Contig * init_contig(char * name, size_t size);

/** Recursively free all memory.
 *
 * This functions calls free_block on each Block in its block field.
 *
 * If the Contig has an IntervalTree defined, it will free it with free_interval_tree.
 *
 * @param contig pointer to a contig, may be NULL
 * */
void free_contig(Contig * contig);

/** Recursively print contig. */
void print_contig(Contig * contig);

/** Find index of downstream Block nearest the query point */
uint anchor(Contig * contig, uint x);

/** Given two points, find all blocks overlapping them */
Contig * get_overlapping(Contig * contig, uint a, uint b);

/** Given two points, find the number of blocks they overlap */
uint count_overlaps(Contig * contig, uint a, uint b);

Contig * get_flanks(Contig * contig, size_t n, bool left);

/** Sort Block objects by start position */
void sort_blocks_by_start(Contig * contig);

/** Sort Block objects by stop position */
void sort_blocks_by_stop(Contig * contig);

#endif
