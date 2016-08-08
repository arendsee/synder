#ifndef __CONTIGUOUS_H__
#define __CONTIGUOUS_H__

#include <stdlib.h>
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

typedef enum {
  SI_GOOD  = 0, // Determinable boundries (case A,B,C,D)
  SI_LSCR  = 1, // Left hand undeterminable (case F)
  SI_RSCR  = 2, // Right hand undeterminable (case F)
  SI_BSCR  = 3, // Both undeterminable (case F on both sides)
  SI_LBND  = 4, // Left unbound, case E to the left of nearest block
  SI_RBND  = 5, // Right unbound, case E to the right of nearest block
  SI_LMOST = 6, // Query is further left than any syntenic block
  SI_RMOST = 7, // Query is further right than any syntenic block
  SI_BMOST = 8, // Query extends beyond all syntenic blocks on both ends
  SI_ERROR = -1 // Something is extremely wrong
} SearchIntervalFlag;

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
void contiguous_query(Synmap * syn, FILE * intfile, bool pblock);

#endif
