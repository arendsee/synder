#ifndef __BLOCK_H__
#define __BLOCK_H__

#include "global.hpp"

Block *init_Block(long, long);

/** A convenience function for setting some of the variables in Block */
void set_Block(
    Block*  block,
    long    start,
    long    stop,
    double  score,
    char    strand,
    Contig* parent,
    Block*  over,
    size_t  linkid
);

void free_Block(Block *);

/** A clean TAB-delimited output suitable for giving to the user */
void print_Block(Block *);

/** A diagnostic function designed for use in synmap dumps
 *
 * This function prints the homologous pair of blocks, not just the the input
 * block.
 *
 */
void print_verbose_Block(Block * blk);

/** Determine whether interval (a,b) overlaps interval (c,d)
 *
 * @param a1 start of first interval
 * @param a2 stop of first interval
 * @param b1 start of second interval
 * @param b2 stop of second interval
 *
 * @return TRUE if the intervals overlap
 */
bool overlap(long a1, long a2, long b1, long b2);

/** Determine whether two Blocks overlap 
 *
 * @return TRUE if they overlap 
 */
bool block_overlap(Block * a, Block * b);

/** Calculate the length of the overlap of two intervals */
long overlap_length_ll(long a1, long a2, long b1, long b2);

/** Calculate the length of the overlap of two intervals */
long overlap_length(Block * a, Block * b);

/** Remove a block and redirect all links to it
 *
 * WARNING: I have not tested this function
 */
void delete_Block(Block* block);

/** Transform b into (a U b), relink as needed
 *
 * After calling this function, a is removed from the datastructure and
 * connections to it are redirected.
 *
 * WARNING: this function only relinks Block->cor connections
 *
 * @param a Block to be merged (and deleted)
 * @param b Block to hold the final union of a and b
 */
void merge_block_a_into_b(Block * a, Block * b);

/** compare intervals by stop */
int block_cmp_stop(const void *, const void *);

/** compare intervals by start */
int block_cmp_start(const void *, const void *);

#endif
