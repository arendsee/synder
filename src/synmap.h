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

/** Link blocks by next and prev stop and next and prev start */
void link_block_corners(Synmap * syn);

/** Link Contig to first and last blocks*/
void set_contig_corners(Synmap * syn);

/** Set a unique index for each set of overlapping sequences */
void set_overlap_group(Synmap * syn);

/** Link each node to its adjacent neighbors */
void link_adjacent_blocks(Synmap * syn);

/** Link each node to its contiguous neighbor */
void link_contiguous_blocks(Synmap * syn);

/** Checks invariants - dies if anything goes wrong */
void validate_synmap(Synmap * syn);

#endif
