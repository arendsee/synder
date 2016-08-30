#ifndef __MAP_H__
#define __MAP_H__

#include <stdint.h>

#include "global.h"
#include "block.h"
#include "synmap.h"

#define ANCHORED  0
#define BOUND     1
#define UNBOUND   2
#define EXTREME   3
#define BEYOND    4

#define cmloc(cmap,blk) cmap->map[block->linkid]

/**
 * @brief Print target regions from a given query
 *
 * Generates a new contiguous map, then searches map for the target regions as
 * passed in through the intfile. 
 *
 * @param Synmap* syn Synteny db
 * @param FILE* intfile GFF file containing query-side target regions, may be
 *                      streamed from STDIN
 * @param bool pblock Toggle to print flanking query and target side blocks to
 *                    a given search region
 * 
 */
void find_search_intervals(Synmap * syn, FILE * intfile, bool pblock);

double calculate_score(long start, long stop, Block * blk);

double calculate_target_score(long a1, long a2, Block * bounds[2]);

#endif
