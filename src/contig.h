#ifndef __CONTIG_H__
#define __CONTIG_H__

#include "global.h"
#include "block.h"
#include "contiguous_set.h"
#include "interval_tree.h"
#include "itree_result.hpp"

#include <array>
#include <vector>
#include <string>
#include <list> 

#include "global.h"
#include "block.h"
#include "synmap.h"
#include "contig.h"
#include "contig_result.h"

class Contig {
private:
    const static long default_length = 1000000000;

    double calculate_score(long start, long stop, Block* blk);

    double calculate_target_score(long a1, long a2, Block* bounds[2]);

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

    /**
     * @brief Print target regions from a given query
     *
     * Generates a new contiguous map, then searches map for the target regions as
     * passed in through the intfile. 
     *
     * @param Synmap* syn Synteny db
     * @param FILE* intfile GFF file containing query-side target regions, may be
     *                      streamed from STDIN
     * 
     */
    void find_search_intervals(Bound& bounds, char* seqname);
};

#endif
