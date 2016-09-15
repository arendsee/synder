#ifndef __SYNMAP_H__
#define __SYNMAP_H__

#include "global.h"
#include "genome.h"
#include "global.h"
#include "contiguous_set.h"
#include "arguments.h"

#include <iterator>
#include <list>


/** A pair of syntenically linked Genome objects  */
class Synmap
{
private:

    Genome* genome[2];
    FILE* synfile;
    FILE* tclfile;
    FILE* qclfile;
    int swap;
    long k;
    char trans;

    // loads synfile and calls the below functions in proper order
    void load_blocks();

    // wrappers for Genome functions
    void set_contig_lengths();
    void link_block_corners();
    void set_contig_corners();
    void merge_doubly_overlapping_blocks();
    void set_overlap_group();
    void link_adjacent_blocks();
    void link_contiguous_blocks(long k);

    /** Checks invariants - dies if anything goes wrong */
    void validate();

public:

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
    Synmap(Arguments& args);
    ~Synmap();

    Contig* get_contig(size_t gid, char* contig_name);

    /** Recursively print a synteny map. */
    void print(bool forward=true);

    void dump_blocks();

    /** Count blocks overlapping intervals in intfile
     *
     * Prints the input sequence name and count to STDOUT in TAB-delimited format.
     *
     * @param syn synmap, where the query and gff_file reference the same genome.
     * @param gff_file GFF format file, 9th column is treated as the interval name.
     */
    void count(FILE* gff_file);

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
    void map(FILE* gff_file);

    void find_search_intervals(FILE* intfile);

};

#endif
