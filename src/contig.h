#ifndef __CONTIG_H__
#define __CONTIG_H__

#include "global.h"
#include "block.h"
#include "itree/itree.h"
#include "itree/search.h"

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
Contig *init_Contig(char *name, size_t size);

/** Recursively free all memory.
 *
 * This functions calls free_Block on each Block in its block field.
 *
 * If the Contig has an IntervalTree defined, it will free it with free_IntervalTree.
 *
 * @param contig pointer to a contig, may be NULL
 * */
void free_Contig(Contig * contig);

/** Recursively free all memory EXCEPT blocks
 *
 * Like free_Contig but does not free all the Blocks. This is useful for
 * structures that hold blocks which do not belong to them (i.e. pointers to
 * blocks held elsewhere).
 *
 * @param contig pointer to a contig, may be NULL
 * */
void free_partial_Contig(Contig * contig);

/** Recursively print contig. */
void print_Contig(Contig * contig, bool forward);

/** A wrapper for Contig that includes IntervalResult flags
 *
 *  From ResultInterval:
 *  1. inbetween - query is between (not overlapping) two search intervals
 *  2. leftmost - query is further left than any interval
 *  3. rightmost - query is further right than any interval
 */
typedef struct ResultContig{
    Contig * contig;
    bool inbetween;
    bool leftmost;
    bool rightmost;
} ResultContig;

ResultContig * init_ResultContig(Contig *, IntervalResult *);

void free_ResultContig(ResultContig *);

void free_partial_ResultContig(ResultContig *);

void print_ResultContig(ResultContig *);

/** Find index of downstream Block nearest the query point */
uint anchor(Contig * contig, uint x);

/** Given two points, find all blocks overlapping or flanking them
 *
 * If there is at least one overlapping block, return all overlapping blocks
 *
 * Otherwise, return the blocks above and below the input region
 *
 * If there is only one flanking interval (i.e., the query is beyond any
 * syntenic interval), return just the nearest interval.
 *
 * */
ResultContig *get_region(Contig * contig, uint a, uint b);

/** Given two points, find the number of blocks they overlap */
uint count_overlaps(Contig * contig, uint a, uint b);

// /** Get intervals flanking input interval */
// Contig * get_flanks(Contig * contig, size_t n, bool left);

/** Sort Block objects by start position */
void sort_blocks_by_start(Contig * contig);

/** Sort Block objects by stop position */
void sort_blocks_by_stop(Contig * contig);

/** Find closest block above/below a given value */
Block * closest_block(Contig * con, int x, Direction d);

#endif
