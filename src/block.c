#include <stdio.h>
#include <stdlib.h>

#include "block.h"

/** Initialize values in a block */
void set_Block(
    Block* block,
    long start,
    long stop,
    float score,
    char strand,
    Contig* parent,
    Block* over)
{
    block->pos[0] = start;
    block->pos[1] = stop;
    block->score = score;
    block->strand = strand;
    block->parent = parent;
    block->over = over;
}

/** Allocate memory for a block and set each field.
 *
 * @param start  query start position
 * @param stop   query stop position
 *
 * @return pointer to new Block
 *
 */
Block* init_Block(long start, long stop)
{
    Block* block = (Block*)calloc(1, sizeof(Block));
    block->pos[0] = start;
    block->pos[1] = stop;
    block->strand = '.';
    return (block);
}

/** Remove this block from the web of linkage
 *
 * The block is not freed, since it is assumed the memory is in a memory pool
 * managed by a Contig object.
 */
void delete_Block_directed(Block* block, Direction d)
{
    // TODO - double check this logic

    // Transformed indices for Block->cor
    // ----------------------------------
    // a   b          c   d
    // <--->   <--->  <--->
    //   a - previous element by start
    //   b - previous element by stop
    //   c - next element by start
    //   d - next element by stop
    int idx_a = (!d * 2) + !d; // - 0
    int idx_b = ( d * 2) + !d; // - 2
    int idx_c = (!d * 2) +  d; // - 1
    int idx_d = ( d * 2) +  d; // - 3

    if (block == block->parent->cor[d])
        block->parent->cor[d] = block->cor[idx_d];

    if (block == block->parent->cor[2 + d])
        block->parent->cor[2 + d] = block->cor[2 + idx_b];

    if (block->cor[idx_a] != NULL)
        block->cor[idx_a]->cor[idx_c] = block->cor[idx_c];

    if (block->cor[idx_b] != NULL)
        block->cor[idx_b]->cor[idx_d] = block->cor[idx_d];

    if (block->cnr[d] != NULL)
        block->cnr[d]->cnr[!d] = block->cnr[!d];

    if (block->adj[d] != NULL) {
        // Normally we link block->adj[1]->adj[0] to block->cor[PREV_STOP], but
        // in the case below, we need to link to NEXT_STOP instead:
        //                    <====>   block->adj[1]
        // i+1       <=====>           block->cor[NEXT_STOP]
        // i         <=====>           block
        // i-1     <====>              block->cor[PREV_STOP]
        if (block->cor[idx_d]->pos[d] == block->cor[idx_c]->pos[d]) {
            block->adj[d]->adj[!d] = block->cor[idx_d];
        } else {
            block->adj[d]->adj[!d] = block->cor[idx_c];
        }
    }
}
void delete_Block(Block* block)
{
    delete_Block_directed(block, LO);
    delete_Block_directed(block, HI);

    // set this to prevent delection of homolog from looping
    block->parent = NULL;
    if (block->over->parent != NULL) {
        delete_Block(block->over);
    }
}

/** Print all fields in this block (TAB-delimited). */
void print_Block(Block* block)
{
    printf("%s\t%zu\t%zu\t%s\t%zu\t%zu\t%c\n",
        block->parent->name,
        block->pos[0] + global_out_start,
        block->pos[1] + global_out_stop,
        block->over->parent->name,
        block->over->pos[0] + global_out_start,
        block->over->pos[1] + global_out_stop,
        block->strand);
}

/**
 * Determine whether interval (a,b) overlaps interval (c,d)
 *
 * @param a1 start of first interval
 * @param a2 stop of first interval
 * @param b1 start of second interval
 * @param b2 stop of second interval
 *
 * @return TRUE if the intervals overlap
 */
bool overlap(long a1, long a2, long b1, long b2)
{
    return a1 <= b2 && a2 >= b1;
}

/**
 * Determine whether two Blocks overlap 
 *
 * @return TRUE if they overlap 
 */
bool block_overlap(Block* a, Block* b)
{
    return a->pos[0] <= b->pos[1] && a->pos[1] >= b->pos[0];
}

/** Compare by Block stop position */
int block_cmp_stop(const void* ap, const void* bp)
{
    Block* a = *(Block**)ap;
    Block* b = *(Block**)bp;
    return (int)(a->pos[1] > b->pos[1]) - (int)(a->pos[1] < b->pos[1]);
}

/** Compare by Block start position */
int block_cmp_start(const void* ap, const void* bp)
{
    Block* a = *(Block**)ap;
    Block* b = *(Block**)bp;
    return (int)(a->pos[0] > b->pos[0]) - (int)(a->pos[0] < b->pos[0]);
}

long get_set_bound(Block* blk, Direction d)
{
    if (blk->cnr[d] != NULL) {
        return get_set_bound(blk->cnr[d], d);
    }
    return blk->pos[d];
}
