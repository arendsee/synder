#ifndef __SYNMAP_ANALYSIS_H__
#define __SYNMAP_ANALYSIS_H__

#include <stdio.h>

#include "synmap.h"

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
