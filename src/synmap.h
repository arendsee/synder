#ifndef __SYNMAP_H__
#define __SYNMAP_H__

#include "global.h"
#include "genome.h"
#include "global.h"
#include "contiguous_set.h"

#include <iterator>
#include <list>


/** A pair of syntenically linked Genome objects  */
class Synmap
{
private:
    void check_args(size_t line_no, size_t nargs, size_t correct_nargs);

    /** Link blocks by next and prev stop and next and prev start */
    void link_block_corners();

    /** Link Contig to first and last blocks */
    void set_contig_corners();

    /** Set a unique index for each set of overlapping sequences */
    void set_overlap_group();

    /** Merge blocks that overlap on both query and target sides
     *
     * For now, I will use a weighted average for the merged Block score.
     *
     */
    void merge_all_doubly_overlapping_blocks();

    /** Link each node to its adjacent neighbors
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
    void link_adjacent_blocks();
    /** Do one half of the job */
    void link_adjacent_blocks_directed(Contig* con, Direction d);

    /** Link each node to its contiguous neighbor */
    void link_contiguous_blocks(long k);

    /** Build and link ContiguousSet structures */
    void build_contiguous_sets();

public:

    Genome * genome[2];

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
    Synmap(FILE * synfile, int swap, long k, char trans, bool validate_on);

    ~Synmap();

    // syn->genome[(gid)]->contig[(cid)]
    Contig * get_contig(size_t gid, size_t cid);

    // syn->genome[(gid)]
    Genome * get_genome(size_t gid);

    /** Recursively print a synteny map. */
    void print(bool forward);

    void dump_blocks();

    /** Checks invariants - dies if anything goes wrong */
    void validate();

    /** Count blocks overlapping intervals in intfile
     *
     * Prints the input sequence name and count to STDOUT in TAB-delimited format.
     *
     * @param syn synmap, where the query and gff_file reference the same genome.
     * @param gff_file GFF format file, 9th column is treated as the interval name.
     */
    void count(FILE * gff_file);

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
    void map(FILE * gff_file);

};

#endif
