#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "contig.h"

Contig *init_Contig(char *name, size_t size, long length)
{
  Contig *con  = (Contig *) malloc(sizeof(Contig));
  con->name    = strdup(name);
  con->length  = length;
  con->size    = size;
  con->block   = (Block *) calloc(size, sizeof(Block));
  con->head[0] = NULL;
  con->head[1] = NULL;
  con->tail[0] = NULL;
  con->tail[1] = NULL;
  con->itree   = NULL;
  return con;
}

void free_Contig(Contig * contig)
{
  if (contig != NULL) {
    if(contig->block){
        free(contig->block);
    }
    if (contig->itree != NULL){
      free_IntervalTree(contig->itree);
    }
    free(contig->name);
    free(contig);
  }
}

void print_Contig(Contig * contig, bool forward)
{
  printf("%lu\t%s\n", contig->size, contig->name);
  Block * blk = contig->head[!forward];
  Corner d = forward ? NEXT_START : NEXT_STOP;
  for(; blk != NULL; blk = blk->cor[d]){
    print_Block(blk);
  }
}

ResultContig * init_ResultContig(Contig * contig, IntervalResult * ir)
{
  ResultContig * rc = (ResultContig *)malloc(sizeof(ResultContig));
  rc->size          = ir->iv->size;
  rc->block         = (Block **)malloc(rc->size * sizeof(Block*));
  rc->contig        = contig;
  rc->inbetween     = ir->inbetween;
  rc->leftmost      = ir->leftmost;
  rc->rightmost     = ir->rightmost;
  for (size_t i = 0; i < rc->size; i++) {
    rc->block[i] = (Block*)ir->iv->v[i].link;
  }
  return rc;
}

void free_ResultContig(ResultContig * rc)
{
  if(rc != NULL){
    if(rc->block != NULL){
        free(rc->block);
    }
    free(rc);
  }
}

/*
void print_ResultContig(ResultContig * rc)
{
    printf("inbetween=%i leftmost=%i rightmost=%i\n", rc->inbetween, rc->leftmost, rc->rightmost);
    print_Contig(rc->contig, true);
}
*/

// Map blocks to the Interval structures used in itree
IA *ia_from_blocks(Contig * con)
{
  Block * blk = con->head[0];

  // Count the number of blks in the linked list
  // Since blocks can be deleted, we cannot just use con->size
  size_t n = 0;
  for(; blk != NULL; blk = blk->cor[1]){ n++; }

  IA *ia = init_set_IA(n);
  // Array index for ia->v array of intervals
  size_t i = 0;
  // Reset blk, since was wound to NULL above
  blk = con->head[0];
  // Map the Block list to the new IA structure
  for (; blk != NULL; blk = blk->cor[1], i++) {
    ia->v[i].start = blk->pos[0];
    ia->v[i].stop  = blk->pos[1];
    ia->v[i].link  = (void *) blk;
  }
  return ia;
}

ResultContig *get_region(Contig * con, long a, long b)
{
  // Build itree if necessary
  if (con->itree == NULL)
    con->itree = build_tree(ia_from_blocks(con));

  // Search itree
  Interval * inv = init_Interval(a, b);
  IntervalResult *res = get_interval_overlaps(inv, con->itree);
  free(inv);

  // itree returns the flanks for queries that overlap nothing. However, I
  // need all the intervals that overlap these flanks as well.
  IntervalResult *tmp_a = NULL;
  IntervalResult *tmp_b = NULL;
  if(res->inbetween){
      // If inbetween, itree should have returned the two flanking blocks
      if(res->iv->size == 2){
          tmp_a = get_interval_overlaps(&res->iv->v[0], con->itree);
          tmp_b = get_interval_overlaps(&res->iv->v[1], con->itree);
          free_IV(res->iv);
          res->iv = tmp_a->iv;
          merge_IV(tmp_a->iv, tmp_b->iv);
          free(tmp_a);
          free_IntervalResult(tmp_b);
      } else {
          fprintf(stderr, "itree is broken, should return exactly 2 intervals for inbetween cases\n");
          exit(EXIT_FAILURE);
      }
  }
  else if(res->leftmost || res->rightmost){
      if(res->iv->size == 1){
          tmp_a = get_interval_overlaps(&res->iv->v[0], con->itree);
          free_IV(res->iv);
          res->iv = tmp_a->iv;
          free(tmp_a);
      } else {
          fprintf(stderr, "itree is broken, should return only 1 interval for left/rightmost cases\n");
          exit(EXIT_FAILURE);
      }
  }

  ResultContig * resultcontig = init_ResultContig(con, res);

  free_IntervalResult(res);

  return resultcontig;
}

long count_overlaps(Contig * con, long a, long b)
{
  if (con->itree == NULL)
    con->itree = build_tree(ia_from_blocks(con));
  Interval *inv = init_Interval(a, b);
  long count = count_interval_overlaps(inv, con->itree);
  free(inv);
  return count;
}

void sort_blocks(Contig * contig, bool by_stop)
{
  if (contig != NULL && contig->block != NULL) {
    if(by_stop) {
        qsort(contig->block, contig->size, sizeof(Block *), block_cmp_stop);
    } else {
        qsort(contig->block, contig->size, sizeof(Block *), block_cmp_start);
    }
  } else {
    fprintf(stderr, "Contig or contig block not defined\n");
  }
}
