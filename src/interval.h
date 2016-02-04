#include "ftypes.h"

/**
 * This structure currently only contains the most basic return values.
 * Later I my expand this to cover whatever other metrics I gather.
 */
typedef struct {
    uint qseqid;
    uint qstart;
    uint qstop;
    uint tseqid;
    uint tstart;
    uint tstop;
    bool mismatch;
} MapResult;

/**
 * Determine wither interval (a,b) overlaps interval (c,d)
 */
bool overlap(uint, uint, uint, uint);

/**
 * Find index of downstream Block nearest the query point
 */
uint anchor(uint, Contig *);

/**
 * Given two points, find all blocks overlapping them
 */
Contig * get_overlapping(uint, uint, Contig *);

/**
 * Given two points, find the blocks flanking them
 */
Contig * get_flanks(uint, uint, Contig *, uint, uint);

/**
 * Given two points, find the number of blocks they overlap
 */
uint count_overlaps(uint, uint, Contig *);

/**
 * Given two points, get the expected location on the target
 */
MapResult map(uint, uint, Contig *);
