#ifndef __CONTIGUOUS_SET_H__
#define __CONTIGUOUS_SET_H__

#include "global.h"
#include "interval.hpp"
#include "feature.h"
#include "linked_interval.hpp"
#include "interval_tree.hpp"
#include "block.h"

#include <array>
#include <algorithm>


/** Contiguous set of non-overlapping adjacent homologous pairs of Blocks */
class ContiguousSet : public LinkedInterval<ContiguousSet>, public Interval<ContiguousSet>
{
private:
    void add_side_(Block* blk);
    bool blocks_conflict(Block* a, Block* b);
    void ContiguousSet_side_(Block* blk);

public:

    size_t               size = 0;
    // TODO Use LinkedInterval cor
    ContiguousSet*       next = nullptr;
    ContiguousSet*       prev = nullptr;
    std::array<Block*,2> ends = {{ nullptr }};

    ContiguousSet();

    /** Homolog constructor */
    ContiguousSet(ContiguousSet* cset);
    ContiguousSet(Block* blk, size_t t_id);

    Block* front() { return ends[0]; }
    Block* back()  { return ends[1]; }

    /** Destructor for a single ContiguousSet in a linked list
     *
     * This function preserves the ContiguousSet linked list
     *  - contig->end are reset, if needed
     *  - ContiguousSet->next and ContiguousSet->prev are relinked
     *  - ContiguousSet->over is freed as well
     *  - Free the parent ContiguousSet IntervalTree, since it will become corrupt
     *
     */
    ~ContiguousSet();

    /** Determine whether two blocks are contiguous
     *
     * Each block consists of query and target intervals
     *
     * Let these 4 resulting intervals be aq, at, bq and bt
     *
     * Each interval is defined by 3 variables:
     *   1. s - target strand
     *   2. c - target chromosome/scaffold name
     *   3. g - non-overlapping group id
     *
     * I will identify each of these variables by appending [scg] to the
     * interval name, e.g. ats or bqg.
     *
     * The blocks are contiguous if and only if all of the following are true
     *   1. aqs == bqs
     *   2. ats == bts
     *   3. aqc == bqc
     *   4. atc == btc
     *
     *   #1 will always be true, since strand is relative to query.
     */
    bool are_contiguous(Block* blk_a, Block* blk_b, long k);

    /** Add a new block to the set if it is contiguous
     *
     * @return bool - true if block is contiguous and was added, false otherwise
     */
    bool add_block(Block* blk, long k);

    /** Add a block without checking for contiguity */
    void force_add_block(Block* blk);

    static bool strictly_forbidden(Block* a, Block* b, long k);
};

#endif
