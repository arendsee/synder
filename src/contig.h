#ifndef __CONTIG_H__
#define __CONTIG_H__

#include <stdlib.h>
#include <stdbool.h>

#include "block.h"

#include "itree/itree.h"
#include "itree/search.h"

#ifndef uint
#define uint unsigned int
#endif

#define NEXT_BLOCK_BYSTART(con, blk) (blk->startid + 1) < con->size ? con->block[blk->startid + 1] : NULL

#define PREV_BLOCK_BYSTOP(con, blk) (blk->stopid > 0) ? con->by_stop[blk->stopid - 1] : NULL

#define GET_RESULT_BLOCK(result, i) (Block*)result->iv->v[(i)].link

#define CB(con, i) (con)->block[(i)]

#define CB_STOP(con, i)  (con)->block[(i)]->stop
#define CB_START(con, i) (con)->block[(i)]->start

#define CB_STOPID(con, i)  (con)->block[(i)]->stopid
#define CB_STARTID(con, i) (con)->block[(i)]->startid

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
  char *name;
  struct IntervalTree *itree;
  size_t size;
  Block **block;
  Block **by_stop;
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

void print_ResultContig(ResultContig *);

/** Find index of downstream Block nearest the query point */
uint anchor(Contig * contig, uint x);

/** Given two points, find all blocks overlapping or flanking them
 *
 * If there is at least one overlapping block, return all overlapping blocks
 *
 * Otherwise, return the blocks above and below the input region
 *
 * If there is no block above or below, return NULL
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

#endif
