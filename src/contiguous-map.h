#ifndef __CONTIGUOUS_MAP_H__
#define __CONTIGUOUS_MAP_H__

#include <stdlib.h>
#include <stdint.h>

#include "global.h"
#include "block.h"
#include "synmap.h"

/** 
 * @brief Contiguous list object holding current Block and its neighbors
 * 
 *
 * Fields:
 * feature  - Current object (query block) 
 * match    - Target match of current feature
 * adj      - Adjacent ContiguousNodes
 * qblkid   - Id of feature (index in parent Contig->block** array)
 * setid    - Unique id for the node's contiguous set
 *
 */
typedef struct ContiguousNode {
  Block *feature;
  Block *match;
  struct ContiguousNode * adj[2];
  size_t qblkid;
  int setid;
} ContiguousNode;

/**
 * @ Initialize a ContiguousNode struct
 */
ContiguousNode * init_ContiguousNode(Synmap * syn, size_t conid, size_t blkid);

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
 * @brief Get smallest or largest value
 *
 * @param cnode
 * @param direction 0 or 1, for finding smallest and largest values, respectively
 */
uint get_set_bound(ContiguousNode * cnode, Direction direction);

/**
 * @brief Debug printing function
 *
 */
void print_ContiguousMap(ContiguousMap * cmap);


/**
 * @brief Debug printing function
 *
 */
void print_ContiguousNode(ContiguousNode * cnode);


/**
 * @brief Populate a new contiguous map from a synteny db
 *
 * @param Synmap* syn Synteny db
 *
 * @return ContiguousMap* populated contiguous map
 */
ContiguousMap *populate_contiguous_map(Synmap * syn);

#endif
