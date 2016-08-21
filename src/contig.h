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
 * @param total sequence length
 *
 * @return pointer to a new Contig
 *
 */
Contig *init_Contig(char *name, size_t size, long length);

/** Recursively free all memory.
 *
 * This functions calls free_Block on each Block in its block field.
 *
 * If the Contig has an IntervalTree defined, it will free it with free_IntervalTree.
 *
 * @param contig pointer to a contig, may be NULL
 */
void free_Contig(Contig * contig);

/** Recursively free all memory EXCEPT blocks
 *
 * Like free_Contig but does not free all the Blocks. This is useful for
 * structures that hold blocks which do not belong to them (i.e. pointers to
 * blocks held elsewhere).
 *
 * @param contig pointer to a contig, may be NULL
 */
void free_partial_Contig(Contig * contig);

/** Recursively print contig. */
void print_Contig(Contig * contig, bool forward);

/** A wrapper for Contig that includes IntervalResult flags
 *
 * Fields
 *  - contig - a pointer to the original Contig
 *  - size   - number of elements in block
 *  - block  - pointer to selection of Block structs in Contig
 *  - flags
 *    1. inbetween - query is between (not overlapping) two search intervals
 *    2. leftmost  - query is further left than any interval
 *    3. rightmost - query is further right than any interval
 */
typedef struct ResultContig{
    Contig * contig;
    size_t size;
    Block ** block;
    bool inbetween;
    bool leftmost;
    bool rightmost;
} ResultContig;

ResultContig * init_ResultContig(Contig *, IntervalResult *);

void free_ResultContig(ResultContig *);

void print_ResultContig(ResultContig *);

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
ResultContig *get_region(Contig * contig, long a, long b);

/** Given two points, find the number of blocks they overlap */
long count_overlaps(Contig * contig, long a, long b);

// /** Get intervals flanking input interval */
// Contig * get_flanks(Contig * contig, long n, bool left);

/** Sort Block objects */
void sort_blocks(Block ** block, size_t size, bool by_stop);

/** Find closest block above/below a given value */
Block * closest_block(Contig * con, long x, Direction d);

#endif
