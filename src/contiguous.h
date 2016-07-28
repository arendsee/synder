#ifndef __CONTIGUOUS_H__
#define __CONTIGUOUS_H__

#include <stdlib.h>
#include <stdint.h>
#include "block.h"
#include "synmap.h"

#define cmloc(cmap,blk) cmap->map[block->linkid]

//Build contiguous set as adjacency list attached to hash map

/** 
 * @brief Contiguous list object holding current Block and its neighbors
 * 
 *
 * Fields:
 * feature 	- Current object
 * match 	- Target match of current feature
 * prev  	- Previous contiguous node, if any
 * next		- Next contiguous node, if any 
 * flag		- Used to indicate several state
 * 			 -3 - Part of a transposed block
 * 			 -2 - Left of last contiguous block
 * 			 -1 - Right of last contiguous block	
 * 			  0 - Normal Interval
 * 			  2 - Query Interval overlaps
 * 			  3 - Target Interval overlaps
 *
 */
typedef struct ContiguousNode {
  Block *feature;
  Block *match;
  struct ContiguousNode *prev;
  struct ContiguousNode *next;
  int flag;
  uint32_t qblkid;
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
 * @brief Populate  a new contiguous map from a synteny db
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
