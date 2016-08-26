#ifndef __BLOCK_H__
#define __BLOCK_H__

#include "global.h"

Block *init_Block(long, long);

/** A convenience function for setting some of the variables in Block */
void set_Block(
    Block*  block,
    long    start,
    long    stop,
    float   score,
    char    strand,
    Contig* parent,
    Block*  over,
    size_t  linkid
);

void free_Block(Block *);

/** A clean TAB-delimited output suitable for giving to the user */
void print_Block(Block *);

/** A diagnostic function designed for use in synmap dumps */
void print_verbose_Block(Block *);

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
