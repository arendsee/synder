#ifndef __BLOCK_H__ 
#define __BLOCK_H__

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
} Block;

Block * init_block(uint, uint, uint, uint, uint);

void free_block(Block *);

void print_block(Block *);

bool overlap(uint, uint, uint, uint);

#endif
