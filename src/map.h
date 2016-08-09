#ifndef __MAP_H__
#define __MAP_H__

#include <stdint.h>

#include "global.h"
#include "block.h"
#include "synmap.h"

#define ANCHORED   1  // 0000 0001
#define BOUND      4  // 0000 0100
#define UNBOUND   16  // 0001 0000
#define EXTREME   64  // 0100 0000

#define F_AA   ( ANCHORED  << 1 ) | ANCHORED 
#define F_AB   ( BOUND     << 1 ) | ANCHORED 
#define F_AU   ( UNBOUND   << 1 ) | ANCHORED 
#define F_AX   ( EXTREME   << 1 ) | ANCHORED 
#define F_BA   ( ANCHORED  << 1 ) | BOUND    
#define F_BB   ( BOUND     << 1 ) | BOUND    
#define F_BU   ( UNBOUND   << 1 ) | BOUND    
#define F_BX   ( EXTREME   << 1 ) | BOUND    
#define F_UA   ( ANCHORED  << 1 ) | UNBOUND  
#define F_UB   ( BOUND     << 1 ) | UNBOUND  
#define F_UU   ( UNBOUND   << 1 ) | UNBOUND  
#define F_UX   ( EXTREME   << 1 ) | UNBOUND  
#define F_XA   ( ANCHORED  << 1 ) | EXTREME  
#define F_XB   ( BOUND     << 1 ) | EXTREME  
#define F_XU   ( UNBOUND   << 1 ) | EXTREME  
#define F_XX   ( EXTREME   << 1 ) | EXTREME  

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

#endif
