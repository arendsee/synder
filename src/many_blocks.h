#ifndef __MANY_BLOCKS_H__
#define __MANY_BLOCKS_H__

#include "global.h"
#include "linked_interval.h"
#include "block.h"
#include "interval_set.h"

#include <list>

class ManyBlocks : public IntervalSet<Block>
{
protected:
    Block* cor[4];    

public:
    // Base over-rides
    Block* front();
    Block* back();
    bool   empty();
    size_t size();
    void   clear();

    Block* front(size_t i);
    Block* terminus(size_t i);

    void link_block_corners();
    void set_overlap_group();
    void link_adjacent_blocks_directed(Direction d);
    void link_adjacent_blocks();
    void merge_overlaps();

};

#endif
