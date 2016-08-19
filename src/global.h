#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

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

typedef struct Synmap Synmap;
typedef struct Genome Genome;
typedef struct Contig Contig;
typedef struct Block Block;

/** A pair of syntenically linked Genome objects  */
struct Synmap {
  size_t size; // should always be 2
  Genome **genome;
};

/** A named set of Contig objects*/
struct Genome {
  char *name;
  size_t size;
  Contig **contig;
};

/** Contiguous sequence object containing list of Block structures and an
 * interval tree to search them.
 *
 * Fields:
 * - name    - a unique name for this contig (e.g. "Chr1")
 * - itree   - an IntervalTree for log(n)+m search of overlapping intervals
 * - size    - number of blocks in the block array
 * - length  - total number of bases in the chromosome/scaffold
 * - block   - array of pointers to Block objects sorted by start
 * - by_stop - array of pointers to Block objects sorted by stop
 * - base    - 0 for 0-based, 1 for 1-based
 */
struct Contig {
  char *name;
  struct IntervalTree *itree;
  size_t size;
  long length;
  Block **block;
  Block **by_stop;
  int base;
};

/** Query interval with directions to matching target
 *
 * Fields:
 * - pos - start and stop positions of the interval
 * - over - pointer to Block on the other genome
 * - parent - pointer to the Contig containing this Block
 * - adj - nearest non-overlapping blocks (0 for left block, 1 for right block)
 * - cnr - adjacent members in the Block's contiguous set (may be NULL)
 * - setid - the id of this Block's contiguous set
 * - grpid - an id shared between this Block and all Block's it overlaps
 *
 *   NOTE: setid and grpid are both initialized to 0 in init_Block. 0 is
 *   reserved for an UNSET id. If a 0 value ever appears after synmap is
 *   initialized, it implies a serious bug.
 */
struct Block {
  long pos[2];
  Block * over;
  Contig * parent;  
  Block * adj[2];
  Block * cnr[2]; // contiguous neighbor
  float score;
  size_t setid;
  size_t grpid;
  char strand;
};

#endif
