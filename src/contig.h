#ifndef __CONTIG_H__
#define __CONTIG_H__

#include "global.h"
#include "many_blocks.h"
#include "many_contiguous_sets.h"
#include "feature.h"

#include <array>
#include <vector>
#include <string>
#include <list>

// Forward declarations
class Genome;

class Contig {
private:
    Feature feat;
    ManyBlocks block;
    ManyContiguousSets cset;

public:

    Contig();
    Contig(std::string genome_name, std::string contig_name, long length=1000000000);

    void set_length(long length);

    void print(bool forward=true, bool print_blocks=true);

    /** Print target regions from a given query */
    void find_search_intervals(const Feature& feat);

    /** Write blocks overlapping intervals in intfile
     *
     * Prints the following TAB-delimited columns to STDOUT:
     * - query entry name, this should be unique for input interval
     * - target contig name
     * - target start position
     * - target stop position
     */
    void map(const Feature& feat);

    /** Count blocks overlapping intervals in intfile
     *
     * Prints the input sequence name and count to STDOUT in TAB-delimited format.
     */
    void count(const Feature& feat);

};

#endif
