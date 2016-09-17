#ifndef __CONTIG_H__
#define __CONTIG_H__

#include "global.h"
#include "block.h"
#include "contiguous_set.h"
#include "interval_tree.h"
#include "interval_set.h"
#include "itree_result.hpp"
#include "bound.h"


#include <array>
#include <vector>
#include <string>
#include <list>

class Contig {
private:
    const static long default_length = 1000000000;

    void build_block_itree();
    void build_cset_itree();
    void build_tree();

public:
    Genome* parent;
    std::string name;
    long length;
    IntervalSet<Block> block;
    IntervalSet<ContiguousSet> cset;

    Contig();
    Contig(std::string name, Genome* parent);

    ~Contig();

    void print(bool forward=true, bool print_blocks=true);

    void merge_doubly_overlapping_blocks();

    void link_contiguous_blocks(long k, size_t &setid);

    /** Given two points, find the number of blocks they overlap */
    long count_overlaps(long a, long b);

    /** This should be called whenever an operation corrupts an itree */
    void clear_cset_tree();

    template <class T>
    IntervalResult<T>* get_region(Bound& bound, IntervalTree<T>* tree, bool get_flanks);

    /** Print target regions from a given query */
    void find_search_intervals(Bound& bounds, char* seqname);

    /** Write blocks overlapping intervals in intfile
     *
     * Prints the following TAB-delimited columns to STDOUT:
     * - query entry name, this should be unique for input interval
     * - target contig name
     * - target start position
     * - target stop position
     */
    void map(Bound& bound, char* seqname);

    /** Count blocks overlapping intervals in intfile
     *
     * Prints the input sequence name and count to STDOUT in TAB-delimited format.
     */
    void count(Bound& bound, char* seqname);

    // TODO need this? Or should IntervalSet handle?
    void sort_blocks(bool by_stop);

};

#endif
