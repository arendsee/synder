#ifndef __CONTIG_H__
#define __CONTIG_H__

#include "global.h"
#include "block.h"
#include "result.h"

/** Contiguous sequence object containing list of Block structures and an
 * interval tree to search them.*/
typedef struct {
    char * name;
    size_t size;
    struct IntervalTree * itree;
    Block ** block;
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

/**
 * Find index of downstream Block nearest the query point
 */
uint anchor(Contig * contig, uint x);

/**
 * Given two points, find all blocks overlapping them
 */
Contig * get_overlapping(Contig * contig, uint a, uint b);

/**
 * Given two points, find the blocks flanking them
 */
Contig * get_flanks(Contig * contig, uint a, uint b, uint nup, uint ndown);

/**
 * Given two points, find the number of blocks they overlap
 */
uint count_overlaps(Contig * contig, uint a, uint b);

/**
 * Given two points, get the expected location on the target
 */
Result map(Contig * contig, uint a, uint b);

#endif
