#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <stdlib.h>

#ifndef uint
#define uint unsigned int
#endif

#define REL_GT(x, y, d)   ((d) ? (x) >  (y) : (x) <  (y))
#define REL_LT(x, y, d)   ((d) ? (x) <  (y) : (x) >  (y))
#define REL_LE(x, y, d)   ((d) ? (x) <= (y) : (x) >= (y))
#define REL_GE(x, y, d)   ((d) ? (x) >= (y) : (x) <= (y))

#define SGCB(syn, gid, cid, bid) syn->genome[(gid)]->contig[(cid)]->block[(bid)]
#define SGC(syn, gid, cid)       syn->genome[(gid)]->contig[(cid)]
#define SG(syn, gid)             syn->genome[(gid)]

typedef enum direction { HI = 1, LO = 0 } Direction;

typedef enum genome_idx { QUERY = 0, TARGET = 1 } Genome_idx;

typedef struct Synmap Synmap;
typedef struct Genome Genome;
typedef struct Contig Contig;
typedef struct Block Block;

/** A pair of syntenically linked Genome objects  */
struct Synmap {
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
 * - name  - a unique name for this contig (e.g. "Chr1")
 * - size  - number of blocks in the block array
 * - block - array of pointers to Block objects
 * - itree - an IntervalTree that allows log(n)+m search of overlapping intervals
 *
 * */
struct Contig {
  char *name;
  struct IntervalTree *itree;
  size_t size;
  Block **block;
  Block **by_stop;
  uint start_sorted:1;
  uint stop_sorted:1;
};

/** Query interval with directions to matching target*/
struct Block {
  uint pos[2];
  Block * over;
  Contig * parent;  
  size_t linkid;
  char strand;
};

#endif
