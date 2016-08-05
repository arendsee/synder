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
Block *init_Block(uint start, uint stop, uint oseqid, uint oblkid, char strand)
{
  Block *block = (Block *) malloc(sizeof(Block));
  block->start = start;
  block->stop = stop;
  block->startid = 0;
  block->stopid = 0;
  block->oseqid = oseqid;
  block->oblkid = oblkid;
  block->linkid = 0; // linkids will be set by load_synmap
  block->strand = strand;
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
  printf("%u\t%u\t%u\t%u\t%u\t%lu\t%lu\t%c\n",
         block->start,
         block->stop,
         block->oseqid,
         block->oblkid,
         block->linkid,
         block->startid,
         block->stopid,
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
  return a->start <= b->stop && a->stop >= b->start;
}

/** Compare by Block stop position */
int block_cmp_stop(const void *ap, const void *bp)
{
  Block *a = * (Block **) ap;
  Block *b = * (Block **) bp;
  return (int)(a->stop > b->stop) - (int)(a->stop < b->stop);
}

/** Compare by Block start position */
int block_cmp_start(const void *ap, const void *bp)
{
  Block *a = * (Block **) ap;
  Block *b = * (Block **) bp;
  return (int)(a->start > b->start) - (int)(a->start < b->start);
}
