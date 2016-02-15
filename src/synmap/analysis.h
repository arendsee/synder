#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include <stdio.h>

#include "synmap.h"
#include "contig/contig.h"

/** Count blocks overlapping intervals in intfile
 *
 * Prints the input sequence name and count to STDOUT in TAB-delimited format.
 *
 * @param syn synmap, where the query and gff_file reference the same genome.
 * @param gff_file GFF format file, 9th column is treated as the interval name.
 */
void analysis_count(Synmap * syn, FILE * gff_file);

/** Write blocks overlapping intervals in intfile
 *
 * Prints the following TAB-delimited columns to STDOUT:
 * - query entry name, this should be unique for input interval
 * - target contig name
 * - target start position
 * - target stop position
 *
 * @param syn synmap, where the query and gff_file reference the same genome.
 * @param gff_file GFF format file, 9th column is treated as the interval name.
 */
void analysis_map(Synmap * syn, FILE * gff_file);

#endif
