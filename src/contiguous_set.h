#ifndef __CONTIGUOUS_SET_H__
#define __CONTIGUOUS_SET_H__

#include "global.h"
#include "interval_tree.h"
#include "linked_interval.h"

// Forward declarations
class Block;


/** Contiguous set of non-overlapping adjacent homologous pairs of Blocks */
class ContiguousSet : LinkedInterval<ContiguousSet>
{
private:
    static void add_block_side_(Block* blk);
    static bool blocks_conflict(Block* a, Block* b);
    static void ContiguousSet_side_(Block* blk);

public:
    size_t         size;
    ContiguousSet* next;
    ContiguousSet* prev;
    size_t         id;
    Block*         ends[2];

    ContiguousSet();
    ContiguousSet(Block* blk);

    ContiguousSet* make_pair(Block* blk);

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

    void print();

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

    bool strictly_forbidden(Block* a, Block* b, long k);

    bool add_block(Block* blk, long k);
};

#endif
