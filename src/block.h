#ifndef __BLOCK_H__ 
#define __BLOCK_H__

#include "global.h"

typedef struct {
    uint start;
    uint stop;
    uint oseqid;
    uint oblkid;
    uint linkid;
} Block;

Block * init_block(uint, uint, uint, uint, uint);

void free_block(Block *);

void print_block(Block *);

/**
 * Determine wether interval (a,b) overlaps interval (c,d)
 */
bool overlap(uint, uint, uint, uint);

#endif
