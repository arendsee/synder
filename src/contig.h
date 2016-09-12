#ifndef __CONTIG_H__
#define __CONTIG_H__

#include "global.h"
#include "block.h"
#include "contiguous_set.h"
#include "itree.h"
#include "contig_result.h"

#include <array>
#include <vector>
#include <string>

class Contig {
private:
    IntervalTree * itree;
    IntervalTree * ctree;
    void build_block_itree();
    void build_cset_itree();
    const static long default_length = 1000000000;
public:
    long length;
    std::string name;
    Genome* parent;
    std::array<Block*,4> cor;
    std::vector<Block*> block;
    ContiguousSet* cset;

    Contig();
    Contig(std::string name, Genome* parent);
    ~Contig();

    void merge_doubly_overlapping_blocks();

    void print(bool forward=true, bool print_blocks=true);

    /** Given two points, find all blocks overlapping or flanking them
     *
     * If there is at least one overlapping block, return all overlapping blocks
     *
     * Otherwise, return the blocks above and below the input region
     *
     * If there is only one flanking interval (i.e., the query is beyond any
     * syntenic interval), return just the nearest interval.
     *
     * */
    ResultContig *get_region(long a, long b, bool is_cset);

    /** Given two points, find the number of blocks they overlap */
    long count_overlaps(long a, long b);

    /** This should be called whenever an operation corrupts an itree */
    void clear_cset_tree();

    /** Sort Block objects */
    static void sort_blocks(Block** blocks, size_t size, bool by_stop);

};

#endif
