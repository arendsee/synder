#ifndef __SYNMAP_H__
#define __SYNMAP_H__

#include "global.h"
#include "bound.h"
#include "genome.h"
#include "arguments.h"
#include "linked_interval.h"

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

    /** Reads a GFF file and calls the appropriate command on each line */
    bool process_gff(FILE* intfile, Command cmd);

};

#endif
