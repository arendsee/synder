#ifndef __SYNMAP_H__
#define __SYNMAP_H__

#include "global.h"
#include "genome.h"

/**
 * Allocate memory for a new Synmap.
 *
 * This struct holds the pair of genomes that wil be compared.
 *
 * @return pointer to the new Synmap
 */
Synmap *init_Synmap();

/** Recursively free all memory allocated to the synteny map.
 * 
 * Calls free_genome on both its Genome children.
 *
 * @param pointer to Synmap struct
 *
 * */
void free_Synmap(Synmap *);

/** Recursively print a synteny map. */
void print_Synmap(Synmap *, bool forward);

void dump_blocks(Synmap *);

/** Link blocks by next and prev stop and next and prev start */
void link_block_corners(Synmap * syn);

/** Link Contig to first and last blocks
 *
 * Requires previous run of link_block_corners
 *
 */
void set_contig_corners(Synmap * syn);

/** Set a unique index for each set of overlapping sequences
 *
 * Requires previous run of set_contig_corners
 *
 */
void set_overlap_group(Synmap * syn);

/** Merge blocks that overlap on both query and target sides
 *
 * For now, I will use a weighted average for the merged Block score.
 *
 */
void merge_all_doubly_overlapping_blocks(Synmap * syn);

/** Link each node to its adjacent neighbors
 *
 * Requires previous run of set_overlap_group
 *
 * Link blocks to nearest non-overlapping up and downstream blocks
 *
 * For example, given these for blocks:
 *  |---a---|
 *            |--b--|
 *             |----c----|
 *                     |---d---|
 *                               |---e---|
 * a->adj := (NULL, b)
 * b->adj := (a, e)
 * c->adj := (a, e)
 * d->adj := (a, e)
 * e->adj := (d, NULL)
 */
void link_adjacent_blocks(Synmap * syn);

/** Link each node to its contiguous neighbor
 *
 * Requires previous run of link_adjacent_blocks
 *
 */
void link_contiguous_blocks(Synmap * syn, long k);

/** Build and link ContiguousSet structures
 *
 * Requires previous run of link_contiguous_blocks
 *
 */ 
void build_contiguous_sets(Synmap* syn);

/** Checks invariants - dies if anything goes wrong
 *
 * Assumes previous run of link_contiguous_blocks
 *
 */
void validate_synmap(Synmap * syn);

#endif
