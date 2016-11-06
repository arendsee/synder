#ifndef __MANY_BLOCKS_H__
#define __MANY_BLOCKS_H__

#include "global.h"
#include "linked_interval.h"
#include "block.h"
#include "interval_set.h"

#include <list>
#include <array>
#include <Rcpp.h>


class ManyBlocks : public IntervalSet<Block>
{
public:
    std::array<Block*, 4> cor = {{ nullptr }};

    // Base over-rides
    ManyBlocks();
    ~ManyBlocks();

    Block* front();
    Block* back();
    bool   empty();
    size_t size();

    Block* front(size_t i);
    Block* corner(size_t i);

    void set_corner(size_t i, Block* blk);

    // Called by Contig:set_contig_corners
    void link_corners();

    void link_block_corners();
    void set_overlap_group(long& offset);
    void link_adjacent_blocks_directed(Direction d);
    void link_adjacent_blocks();
    void merge_overlaps();
    void refresh();

};

#endif
