#include <stdlib.h>
#include <stdio.h>

#include "block.h"

/** Allocate memory for a block and set each field.
 *
 * @param start  query start position
 * @param stop   query stop position
 * @param oseqid index of target Contig
 * @param oblkid index of matching Block on target Contig
 * @param linkid index of this Block pair in metadata array
 *
 * @return pointer to new Block
 *
 * */
Block *init_Block(uint start, uint stop)
{
  Block *block  = (Block *) malloc(sizeof(Block));
  block->pos[0] = start;
  block->pos[1] = stop;
  block->over   = NULL;
  block->parent = NULL;
  block->adj[0] = NULL;
  block->adj[1] = NULL;
  block->cnr[0] = NULL;
  block->cnr[1] = NULL;
  block->setid  = 0; // reserve 0 for unset
  block->grpid  = 0; // reserve 0 for unset
  block->strand = '.';
  return (block);
}

/** Free memory allocated to this block.
 *
 * @param block pointer to a Block, may be NULL
 * */
void free_Block(Block * block)
{
  if (block != NULL) {
    free(block);
  }
}

/** Print all fields in this block (TAB-delimited). */
void print_Block(Block * block)
{
  printf("%s\t%u\t%u\t%s\t%u\t%u\t%c\n",
         block->parent->name,
         block->pos[0] + global_out_base,
         block->pos[1] + global_out_base,
         block->over->parent->name,
         block->over->pos[0] + global_out_base,
         block->over->pos[1] + global_out_base,
         block->strand
  );
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
bool overlap(uint a1, uint a2, uint b1, uint b2)
{
  return a1 <= b2 && a2 >= b1;
}

/**
 * Determine whether two Blocks overlap 
 *
 * @return TRUE if they overlap 
 */
bool block_overlap(Block * a, Block * b)
{
  return a->pos[0] <= b->pos[1] && a->pos[1] >= b->pos[0];
}

/** Compare by Block stop position */
int block_cmp_stop(const void *ap, const void *bp)
{
  Block *a = * (Block **) ap;
  Block *b = * (Block **) bp;
  return (int)(a->pos[1] > b->pos[1]) - (int)(a->pos[1] < b->pos[1]);
}

/** Compare by Block start position */
int block_cmp_start(const void *ap, const void *bp)
{
  Block *a = * (Block **) ap;
  Block *b = * (Block **) bp;
  return (int)(a->pos[0] > b->pos[0]) - (int)(a->pos[0] < b->pos[0]);
}

uint get_set_bound(Block * blk, Direction d){
    if(blk->cnr[d] != NULL){
        return get_set_bound(blk->cnr[d], d);
    }
    return blk->pos[d];
}
