#ifndef __CONTIGUOUS_H__
#define __CONTIGUOUS_H__

#include <stdlib.h>
#include <stdint.h>
#include "block.h"
#include "synmap.h"

#define ANCHORED_L 1
#define ANCHORED_R 2
#define BOUND_L    4
#define BOUND_R    8
#define UNBOUND_L  16
#define UNBOUND_R  32
#define EXTREME_L  64
#define EXTREME_R 128
#define F_AA    ANCHORED_L | ANCHORED_R
#define F_AB    ANCHORED_L | BOUND_R
#define F_AU    ANCHORED_L | UNBOUND_R
#define F_AX    ANCHORED_L | EXTREME_R
#define F_BA    BOUND_L    | ANCHORED_R
#define F_BB    BOUND_L    | BOUND_R
#define F_BU    BOUND_L    | UNBOUND_R
#define F_BX    BOUND_L    | EXTREME_R
#define F_UA    UNBOUND_L  | ANCHORED_R
#define F_UB    UNBOUND_L  | BOUND_R
#define F_UU    UNBOUND_L  | UNBOUND_R
#define F_UX    UNBOUND_L  | EXTREME_R
#define F_XA    EXTREME_L  | ANCHORED_R
#define F_XB    EXTREME_L  | BOUND_R
#define F_XU    EXTREME_L  | UNBOUND_R
#define F_XX    EXTREME_L  | EXTREME_R

#define cmloc(cmap,blk) cmap->map[block->linkid]

typedef enum {
  SI_GOOD   = 0, // Determinable boundries (case A,B,C,D)
  SI_LUNDET = 1, // Left hand undeterminable (case F)
  SI_RUNDET = 2, // Right hand undeterminable (case F)
  SI_BUNDET = 3, // Both undeterminable (case F on both sides) 
  SI_LUNBND = 4, // Left unbound, case E to the left of nearest block
  SI_RUNBND = 5  // Right unbound, case E to the right of nearest block
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
