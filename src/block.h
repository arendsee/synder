#ifndef __BLOCK_H__
#define __BLOCK_H__

#include <stdlib.h>
#include <stdbool.h>

#ifndef uint
#define uint unsigned int
#endif

/** Query interval with directions to matching target*/
typedef struct {
  uint start;
  uint stop;
  uint oseqid;
  uint oblkid;
  uint linkid;
  size_t startid;
  size_t stopid;
  char strand;
} Block;

Block *init_block(uint, uint, uint, uint, uint, char);

void free_block(Block *);

void print_block(Block *);

bool overlap(uint, uint, uint, uint);

bool block_overlap(Block *, Block *);

/** compare intervals by stop
 *
 * \todo test sort_contig_by_start and sort_contig_by_stop functions
 */
int block_cmp_stop(const void *, const void *);

/** compare intervals by start */
int block_cmp_start(const void *, const void *);

#endif
