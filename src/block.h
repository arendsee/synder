#ifndef __BLOCK_H__
#define __BLOCK_H__

#include "global.h"
#include "feature.h"
#include "interval.h"
#include "linked_interval.h"

#include <array>

class ContiguousSet;

class Block : public LinkedInterval<Block>, public Interval<Block>
{
friend class ManyBlocks;
private:
    static void unlink(Block* blk, int u, int d);
    static void move_b_to_a(Block* a, Block* b, int u, int d);
    static void replace_edge(Block* a, Block* b, int u, int d);
    static void dissolve_edge(Block* blk, int u, int d);
    static void merge_block_a_into_b_edge_(Block* a, Block* b, int i);

public:
    // adjacent block in contiguous set
    std::array<Block*, 2> cnr = {nullptr, nullptr};
    // contiguous set id
    ContiguousSet* cset = nullptr;

    Block();
    Block(
        long     start,
        long     stop,
        double   score,
        char     strand,
        Feature* parent,
        size_t   linkid
    );
    ~Block();

    /** A clean TAB-delimited output suitable for giving to the user */
    void print();

    /** Transform b into (a U b), relink as needed
     *
     * After calling this function, a is removed from the datastructure and
     * connections to it are redirected.
     *
     * WARNING: this function only relinks Block->cor connections
     *
     * @param a Block to be merged (and deleted)
     * @param b Block to hold the final union of a and b
     */
    static void merge_block_a_into_b(Block* a, Block* b);
};

#endif
