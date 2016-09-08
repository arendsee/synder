#ifndef __BLOCK_H__
#define __BLOCK_H__

#include "global.h"
#include "contig.h"


/** Query interval with directions to matching target
 *
 * Fields:
 * - pos - start and stop positions of the interval
 * - over - pointer to Block on the other genome
 * - parent - pointer to the Contig containing this Block
 * - corner - prev and next by start and stop
 * - adj - nearest non-overlapping blocks (0 for left block, 1 for right block)
 * - cnr - adjacent members in the Block's contiguous set (may be NULL)
 * - set - pointer to this Block's ContiguousSet
 * - grpid - an id shared between this Block and all Block's it overlaps
 *
 * next and prev used when you need to just iterate through all the block, they
 * also allow easy deletion of blocks.
 *
 * adj and cnr are designed for specialized cases where symmetry is important.
 *
 *   NOTE: setid and grpid are both initialized to 0 in init_Block. 0 is
 *   reserved for an UNSET id. If a 0 value ever appears after synmap is
 *   initialized, it implies a serious bug.
 */
struct Block
{
    Contig        *parent; // Contig parent
    Block         *over;   // homologous block in other genome
    Block         *cor[4]; // next and prev elements by start and stop
    Block         *adj[2]; // adjacent non-overlapping block
    Block         *cnr[2]; // adjacent block in contiguous set
    ContiguousSet *cset;   // contiguous set id
    long          pos[2];  // start and stop positions
    double        score;   // score provided by synteny program
    size_t        grpid;   // overlapping group id;
    char          strand;  // strand [+-.]
    size_t        linkid;  // a unique block id used mostly for debugging
};

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
