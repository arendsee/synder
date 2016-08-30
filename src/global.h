#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>
#include <math.h>


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

#define SGCB(syn, gid, cid, bid) syn->genome[(gid)]->contig[(cid)]->block[(bid)]
#define SGC(syn, gid, cid)       syn->genome[(gid)]->contig[(cid)]
#define SG(syn, gid)             syn->genome[(gid)]

#define LINE_BUFFER_SIZE 512
#define NAME_BUFFER_SIZE 128

// A value of 0 or 1 with is added to the starts and stops of all printed intervals
int global_in_start;
int global_in_stop;
int global_out_start;
int global_out_stop;

/**
 * If the input is truly 0-based, but the user says it is 1-based, we can get
 * an overflow if we subtract 1 from 0 (and of course our output will be
 * incorrect). This function checks whether all start and stop positions are
 * greater than 0 if 1-based.
 */
void check_in_offset(size_t start, size_t stop);

typedef enum direction { LO = 0, HI = 1 } Direction;

typedef enum genome_idx { QUERY = 0, TARGET = 1 } Genome_idx;

typedef enum corner {
    PREV_START = 0,   // previous element as ordered by start
    NEXT_START = 1,   // next element as ordered by start
    PREV_STOP  = 2,   // previous element as ordered by stop
    NEXT_STOP  = 3    // next element as ordered by stop
} Corner;

typedef struct Synmap Synmap;
typedef struct Genome Genome;
typedef struct Contig Contig;
typedef struct ContiguousSet ContiguousSet;
typedef struct Block Block;

/** A pair of syntenically linked Genome objects  */
struct Synmap {
  size_t size; // should always be 2
  Genome ** genome;
};

/** A named set of Contig objects*/
struct Genome {
  char * name;
  size_t size;
  Contig ** contig;
};

/** Contiguous sequence object containing list of Block structures and an
 * interval tree to search them.
 *
 * Fields:
 * - name    - a unique name for this contig (e.g. "Chr1")
 * - idx     - a unique index between 0 and (Genome->size - 1)
 * - length  - total number of bases in the chromosome/scaffold
 * - cor     - four terminal Blocks 
 *   [0] first by start
 *   [1] last by start
 *   [2] first by stop
 *   [3] last by stop
 * - size    - total number of Blocks in Contig (including deleted ones)
 * - block   - memory pool for all Block structs
 * - itree   - an IntervalTree for log(n)+m search of overlapping intervals (Blocks)
 * - ctree   - an IntervalTree for searching ContiguousSets
 */
struct Contig {
  char * name;
  size_t idx;
  // handle Blocks
  long length;
  Block * cor[4];
  size_t size;
  Block * block;
  struct IntervalTree * itree;
  // handle ContiguousSets
  ContiguousSet * cset;
  struct IntervalTree * ctree;
};

/** Contiguous set of non-overlapping adjacent homologous pairs of Blocks */
struct ContiguousSet {
    long bounds[2];
    size_t size;
    Contig * parent;
    ContiguousSet * next;
    ContiguousSet * prev;
    ContiguousSet * over;
    char strand;
    size_t id; // mostly for debugging
    Block * ends[2];
};

/** Query interval with directions to matching target
 *
 * Fields:
 * - pos - start and stop positions of the interval
 * - over - pointer to Block on the other genome
 * - parent - pointer to the Contig containing this Block
 * - corner - prev and next by start and stop
 * - adj - nearest non-overlapping blocks (0 for left block, 1 for right block)
 * - cnr - adjacent members in the Block's contiguous set (may be NULL)
 * - set - pointer to this Block's ContiguousSet
 * - grpid - an id shared between this Block and all Block's it overlaps
 *
 * next and prev used when you need to just iterate through all the block, they
 * also allow easy deletion of blocks.
 *
 * adj and cnr are designed for specialized cases where symmetry is important.
 *
 *   NOTE: setid and grpid are both initialized to 0 in init_Block. 0 is
 *   reserved for an UNSET id. If a 0 value ever appears after synmap is
 *   initialized, it implies a serious bug.
 */
struct Block {
  Contig        *parent; // Contig parent
  Block         *over;   // homologous block in other genome
  Block         *cor[4]; // next and prev elements by start and stop
  Block         *adj[2]; // adjacent non-overlapping block
  Block         *cnr[2]; // adjacent block in contiguous set
  ContiguousSet *cset;   // contiguous set id
  long          pos[2];  // start and stop positions
  float         score;   // score provided by synteny program
  size_t        grpid;   // overlapping group id;
  char          strand;  // strand [+-.]
  size_t        linkid;  // a unique block id used mostly for debugging
};

#endif
