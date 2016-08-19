#ifndef __BLOCK_H__
#define __BLOCK_H__

#include "global.h"

Block *init_Block(long, long);

void free_Block(Block *);

void print_Block(Block *);

bool overlap(long, long, long, long);

bool block_overlap(Block *, Block *);

/** compare intervals by stop */
int block_cmp_stop(const void *, const void *);

/** compare intervals by start */
int block_cmp_start(const void *, const void *);

/**
 * @brief Get smallest or largest value
 *
 * @param block 
 * @param direction 0 or 1, for finding smallest and largest values, respectively
 */
long get_set_bound(Block * block, Direction direction);

#endif
