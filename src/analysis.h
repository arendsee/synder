#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include <stdio.h>
#include <stdbool.h>

#include "synmap.h"
#include "contig.h"

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

/** Hold a syntenic match independent of a Synmap
 * 
 *  /todo it might be better to actually build the whole input into a Synmap
 *        and then compare Synmap objects.
 *
 */
typedef struct {
    size_t qseqid;
    uint qstart;
    uint qstop;
    size_t tseqid;
    uint tstart;
    uint tstop;
} Link;

/** Identifies syntenic links that agree with the syntenic database
 *
 * @param syn Synmap database
 * @param syn_file FILE pointer to data with synteny format
 */
void analysis_filter(Synmap * syn, FILE * syn_file,
                     bool(*classifier)(Synmap *, Link *, void *),
                     void *);

/** Classify as supported if at least one syn block is nearby
 *
 * This classifier can be passed as the 3rd argument to analysis_filter. It
 * returns true if a block from SYN is within WIDTH characters of the start or
 * start of QUERY in both genomes of SYN.
 *
 * @param syn Synmap database
 * @param query the syntenic block that is being matched against syn
 * @param width a pointer to an unsigned integer describing the distance from
 *              QUERY flanks to search for a parallel block
 */
bool single_advocate(Synmap * syn, Link * query, void * width);

#endif
