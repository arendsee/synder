#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>

#include <cstring>
#include <string>
#include <vector>


#define REL_GT(x, y, d)   ((d) ? (x) >  (y) : (x) <  (y))
#define REL_LT(x, y, d)   ((d) ? (x) <  (y) : (x) >  (y))
#define REL_LE(x, y, d)   ((d) ? (x) <= (y) : (x) >= (y))
#define REL_GE(x, y, d)   ((d) ? (x) >= (y) : (x) <= (y))
#define REL_INC(x, d)     ((d) ? (x++) : (x--))
#define REL_DEC(x, d)     ((d) ? (x--) : (x++))
#define REL_LO_IDX(c, d)  ((d) ? 0 : c->size - 1)
#define REL_HI_IDX(c, d)  ((d) ? c->size - 1 : 0)
#define REL_NEXT(a, i, d) ((d) ? (a)[i+1] : (a)[i-1])
#define REL_ADD(a, b, d)  ((d) ? (a) + (b) : (a) - (b))
#define REL_SUB(a, b, d)  ((d) ? (a) - (b) : (a) + (b))

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

#define LINE_BUFFER_SIZE 512
#define NAME_BUFFER_SIZE 128

// A value of 0 or 1 with is added to the starts and stops of all printed intervals
class Offsets
{
public:
    static int in_start;
    static int in_stop;
    static int out_start;
    static int out_stop;
    // void set_offsets(int offsets[4]){
    //     in_start  = offsets[0];
    //     in_stop   = offsets[1];
    //     out_start = offsets[2];
    //     out_stop  = offsets[3];
    // }
};


/**
 * If the input is truly 0-based, but the user says it is 1-based, we can get
 * an overflow if we subtract 1 from 0 (and of course our output will be
 * incorrect). This function checks whether all start and stop positions are
 * greater than 0 if 1-based.
 */
void check_in_offset(size_t start, size_t stop);

typedef enum direction { LO = 0, HI = 1 } Direction;

typedef enum genome_idx { QUERY = 0, TARGET = 1 } Genome_idx;

typedef enum corner
{
    PREV_START = 0,   // previous element as ordered by start
    NEXT_START = 1,   // next element as ordered by start
    PREV_STOP  = 2,   // previous element as ordered by stop
    NEXT_STOP  = 3    // next element as ordered by stop
} Corner;

#endif
