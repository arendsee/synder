#ifndef __CONTIG_H__
#define __CONTIG_H__

#include "global.h"
#include "block.h"
#include "contiguous_set.h"
#include "itree.h"
#include "itree_result.hpp"

#include <array>
#include <vector>
#include <string>
#include <list> 


class Contig {
private:
    const static long default_length = 1000000000;

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

};

#endif
