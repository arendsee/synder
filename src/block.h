#ifndef __BLOCK_H__
#define __BLOCK_H__

#include "global.h"

Block *init_Block(uint, uint);

void free_Block(Block *);

void print_Block(Block *);

bool overlap(uint, uint, uint, uint);

bool block_overlap(Block *, Block *);

/** compare intervals by stop */
int block_cmp_stop(const void *, const void *);

/** compare intervals by start */
int block_cmp_start(const void *, const void *);

#endif
