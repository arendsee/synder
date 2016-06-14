#ifndef __CONTIGUOUS_H__
#define __CONTIGUOUS_H__

#include <stdlib.h>
#include <stdint.h>
#include "block.h"
#include "synmap.h"

#define cmloc(cmap,blk) cmap->map[block->linkid]

//Build contiguous set as adjacency list attached to hash map

/** Contiguous list object holding current Block and its neighbors
 * 
 *
 * Fields:
 * feature 	- Current object
 * match 	- Target match of current feature
 * prev  	- Previous contiguous node, if any
 * next		- Next contiguous node, if any 
 * flag		- Used to indicate several state
 * 			 -2 - Left of last contiguous block
 * 			 -1 - Right of last contiguous block	
 * 			  0 - Normal Interval
 * 			  2 - Query Interval overlaps
 * 			  3 - Target Interval overlaps
 *
 */


typedef struct  ContiguousNode{
	Block * feature;
	Block * match;
	struct ContiguousNode * prev;
    struct ContiguousNode * next;	
	int flag;
	uint32_t qblkid;
} ContiguousNode;

/** Adjacency list structure to store head and tail 
 * of each contiguous set
 *
 * Fields:
 * first - head of contiguous set
 * last - tail of contiguous set
 * size - number of members of contiguous set
 *
 */
typedef struct ContiguousList {
	ContiguousNode * first;
	ContiguousNode * last;
	size_t  size;
} ContiguousList;

/** Map using linkid to serve as a hashmap to allow quick access
 * of individual nodes of the adjacency list
 *
 * Fields:
 * map - array of pointers to ContiguousNodes
 * size - number of nodes in map
 *
 */
typedef struct ContiguousMap{
	ContiguousNode ** map;
	size_t size;
} ContiguousMap;

/**Create a new contiguous list
 *
 * This initializes an empty contiguous list
 *
 * @return pointer to new contiguous list
 */

ContiguousList * create_contiguous_list();

//Clear existing contiguous list
void free_contiguous_list(ContiguousList *clist);
void free_contiguous_map(ContiguousMap *cmap);

void contiguous_list_push(ContiguousList *clist, ContiguousNode *cnode);
void contiguous_list_reset(ContiguousList *clist, ContiguousNode *cnode);

ContiguousMap * init_contiguous_map(size_t size);
ContiguousMap * populate_contiguous_map(Synmap * syn);

void contiguous_query(Synmap * syn, FILE * intfile);

#endif
