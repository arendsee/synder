#ifndef __SYNMAP_H__
#define __SYNMAP_H__

#include <stdio.h>

#include "genome/genome.h"

#ifndef uint
#define uint unsigned int
#endif

/** A pair of syntenically linked Genome objects  */
typedef struct {
    Genome ** genome;
} Synmap;

/**
 * Allocate memory for a new Synmap.
 *
 * This struct holds the pair of genomes that wil be compared.
 *
 * @return pointer to the new Synmap
 */
Synmap * init_synmap();

/** Recursively free all memory allocated to the synteny map.
 * 
 * Calls free_genome on both its Genome children.
 *
 * @param pointer to Synmap struct
 *
 * */
void free_synmap(Synmap *);

/** Recursively print a synteny map. */
void print_synmap(Synmap *);

/** Build synteny tree from specially formatted file.
 *
 * @warning This function is VERY picky about input. It expects input to be
 * formatted exactly as util/prepare-data.sh produces. You must not feed this
 * function raw synteny files. I currently have no input checks.
 *
 * @param synfile specially formatted synteny file
 *
 * @return pointer to a complete Synmap object
 */
Synmap * load_synmap(FILE *);

/** Count blocks overlapping intervals in intfile
 *
 * Output is printed to STDOUT
 *
 * @param syn synmap, where the query and gff_file reference the same genome.
 * @param gff_file GFF format file, 9th column is treated as the interval name.
 */
void analysis_count(Synmap * syn, FILE * gff_file);

/** Write blocks overlapping intervals in intfile
 *
 * Output is printed to STDOUT
 *
 * @param syn synmap, where the query and gff_file reference the same genome.
 * @param gff_file GFF format file, 9th column is treated as the interval name.
 */
void analysis_map(Synmap * syn, FILE * gff_file);

#endif
