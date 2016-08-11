#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#ifndef uint
#define uint unsigned int
#endif

#define REL_GT(x, y, d)   ((d) ? (x) >  (y) : (x) <  (y))
#define REL_LT(x, y, d)   ((d) ? (x) <  (y) : (x) >  (y))
#define REL_LE(x, y, d)   ((d) ? (x) <= (y) : (x) >= (y))
#define REL_GE(x, y, d)   ((d) ? (x) >= (y) : (x) <= (y))
#define REL_INC(x, d)     ((d) ? (x++) : (x--))
#define REL_DEC(x, d)     ((d) ? (x--) : (x++))
#define REL_LO_IDX(c, d)  ((d) ? 0 : c->size - 1)
#define REL_HI_IDX(c, d)  ((d) ? c->size - 1 : 0)
#define REL_NEXT(a, i, d) ((d) ? (a)[i+1] : (a)[i-1])

#define SGCB(syn, gid, cid, bid) syn->genome[(gid)]->contig[(cid)]->block[(bid)]
#define SGC(syn, gid, cid)       syn->genome[(gid)]->contig[(cid)]
#define SG(syn, gid)             syn->genome[(gid)]

#define LINE_BUFFER_SIZE 512
#define NAME_BUFFER_SIZE 128

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
 * - name    - a unique name for this contig (e.g. "Chr1")
 * - itree   - an IntervalTree for log(n)+m search of overlapping intervals
 * - size    - number of blocks in the block array
 * - length  - total number of bases in the chromosome/scaffold
 * - block   - array of pointers to Block objects sorted by start
 * - by_stop - array of pointers to Block objects sorted by stop
 *
 */
struct Contig {
  char *name;
  struct IntervalTree *itree;
  size_t size;
  unsigned int length;
  Block **block;
  Block **by_stop;
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
 */
struct Block {
  uint pos[2];
  Block * over;
  Contig * parent;  
  Block * adj[2];
  Block * cnr[2]; // contiguous neighbor
  size_t setid;

  // grpid holds the id the overlapping intervals on each genome.  It indexed
  // by linkid. This allows lookup of block adjacency. Two blocks are adjacent
  // if their setids differ by 1 (if on plus strand) or -1 (if on negative
  // strand).
  size_t grpid;
  char strand;
};

#endif
