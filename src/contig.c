#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "contig.h"

Contig *init_Contig(char *name, size_t size, long length)
{
  Contig *con = (Contig *) malloc(sizeof(Contig));
  con->name = strdup(name);
  con->size = size;
  con->length = length;
  con->itree = NULL;
  con->block = (Block **) malloc(size * sizeof(Block *));
  con->by_stop = NULL;
  return con;
}

void free_Contig(Contig * contig)
{
  if (contig != NULL) {
    for (size_t i = 0; i < contig->size; i++) {
      if (contig->block[i] != NULL)
        free_Block(contig->block[i]);
    }
    if (contig->itree != NULL)
      free_IntervalTree(contig->itree);
    if (contig->by_stop != NULL)
      free(contig->by_stop);
    free(contig->block);
    free(contig->name);
    free(contig);
  }
}

void free_partial_Contig(Contig * contig)
{
  if (contig != NULL) {
    if (contig->itree != NULL)
      free_IntervalTree(contig->itree);
    if (contig->by_stop != NULL)
      free(contig->by_stop);
    if(contig->block != NULL)
      free(contig->block);
    free(contig->name);
    free(contig);
  }
}

void print_Contig(Contig * contig, bool forward)
{
  printf("%lu\t%s\n", contig->size, contig->name);
  if(forward){
    for (size_t i = 0; i < contig->size; i++) {
      print_Block(contig->block[i]);
    }
  } else {
    for (size_t i = 0; i < contig->size; i++) {
      print_Block(contig->by_stop[i]);
    }
  }
}

ResultContig * init_ResultContig(Contig * contig, IntervalResult * ir)
{
  ResultContig * rc = (ResultContig *)malloc(sizeof(ResultContig));
  rc->contig        = contig;
  rc->inbetween     = ir->inbetween;
  rc->leftmost      = ir->leftmost;
  rc->rightmost     = ir->rightmost;
  return rc;
}

void free_ResultContig(ResultContig * rc)
{
  if(rc != NULL){
    if(rc->contig != NULL){
        free_Contig(rc->contig);
    }
    free(rc);
  }
}

void free_partial_ResultContig(ResultContig * rc)
{
  if(rc != NULL){
    if(rc->contig != NULL){
        free_partial_Contig(rc->contig);
    }
    free(rc);
  }
}

void print_ResultContig(ResultContig * rc)
{
    printf("inbetween=%i leftmost=%i rightmost=%i\n", rc->inbetween, rc->leftmost, rc->rightmost);
    print_Contig(rc->contig, true);
}

// Map blocks to the Interval structures used in itree
IA *ia_from_blocks(Contig * con)
{
  IA *ia = init_set_IA(con->size);
  for (size_t i = 0; i < con->size; i++) {
    ia->v[i].start = con->block[i]->pos[0];
    ia->v[i].stop = con->block[i]->pos[1];
    ia->v[i].link = (void *) con->block[i];
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

  // Assign returned intervals to Contig
  Contig *contig = init_Contig(con->name, res->iv->size, con->length);
  for (size_t i = 0; i < res->iv->size; i++) {
    contig->block[i] = (Block*)res->iv->v[i].link;
  }

  ResultContig * resultcontig = init_ResultContig(contig, res);

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

void sort_blocks_by_start(Contig * contig)
{
  if (contig != NULL) {
    qsort(contig->block, contig->size, sizeof(Block *), block_cmp_start);
  }
}

void sort_blocks_by_stop(Contig * contig)
{
  if (contig != NULL) {
    if (contig->by_stop == NULL) {
      contig->by_stop = (Block **) malloc(contig->size * sizeof(Block *));
      memcpy(contig->by_stop, contig->block, contig->size * sizeof(Block *));
    }
    qsort(contig->by_stop, contig->size, sizeof(Block *), block_cmp_stop);
  }
}
