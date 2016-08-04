#ifndef __CONTIGUOUS_H__
#define __CONTIGUOUS_H__

#include <stdlib.h>
#include <stdint.h>
#include "block.h"
#include "synmap.h"

#define cmloc(cmap,blk) cmap->map[block->linkid]

// Build contiguous set as adjacency list attached to hash map

typedef enum {
  CNF_UNSET  = -4, // Flag has not yet been set
  CNF_INV    = -3, // Part of a inverted block
  CNF_LEFT   = -2, // Left of last contiguous block
  CNF_RIGHT  = -1, // Right of last contiguous block
  CNF_NORMAL =  0, // Normal Interval
  CNF_QOVER  =  1, // (2) Query Interval overlaps
  CNF_TOVER  =  2, // (3) Target Interval overlaps
  CNF_BOVER  =  3  // Both sides overlaps
} ContiguousNodeFlag;

typedef enum {
  SI_GOOD   = 0, // Determinable boundries (case A,B,C,D)
  SI_LUNDET = 1, // Left hand undeterminable (case F)
  SI_RUNDET = 2, // Right hand undeterminable (case F)
  SI_BUNDET = 3, // Both undeterminable (case F on both sides) 
  SI_LUNBND = 4, // Left unbound, case E to the left of nearest block
  SI_RUNBND = 5  // Right unbound, case E to the right of nearest block
} SearchIntervalFlag;


/** 
 * @brief Contiguous list object holding current Block and its neighbors
 * 
 *
 * Fields:
 * feature - Current object (query block) 
 * match   - Target match of current feature
 * prev    - Previous contiguous node, if any
 * next    - Next contiguous node, if any 
 * qblkid  - Id of feature (index in parent Contig->block** array)
 * setid   - Unique id for the node's contiguous set
 *
 */
typedef struct ContiguousNode {
  Block *feature;
  Block *match;
  struct ContiguousNode *prev;
  struct ContiguousNode *next;
  ContiguousNodeFlag flag;
  size_t qblkid;
  size_t setid;
} ContiguousNode;

/** Map using linkid to serve as a hashmap to allow quick access
 * of individual nodes of the adjacency list
 *
 * Fields:
 * map - array of pointers to ContiguousNodes
 * size - number of nodes in map
 *
 */
typedef struct ContiguousMap {
  ContiguousNode **map;
  size_t size;
} ContiguousMap;


/**
 * @brief Frees contiguous map from memory
 *
 * @param ContiguousMap* cmap The contiguous map to free
 */
void free_ContiguousMap(ContiguousMap * cmap);

/**
 * @brief Initialize a new contiguous map
 *
 * @param size_t size Number of nodes in the contiguous map
 *
 * @return ContiguousMap* New contiguous map
 */
ContiguousMap *init_ContiguousMap(size_t size);

/**
 * @brief Populate a new contiguous map from a synteny db
 *
 * @param Synmap* syn Synteny db
 *
 * @return ContiguousMap* populated contiguous map
 */
ContiguousMap *populate_contiguous_map(Synmap * syn);

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
