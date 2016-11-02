#include "many_blocks.h"

ManyBlocks::ManyBlocks() { }

ManyBlocks::~ManyBlocks(){ }

Block* ManyBlocks::front()
{
    return cor[0];
}

Block* ManyBlocks::corner(size_t i)
{
    if(i > 3) {
        throw "Illegal index for Block::cor in ManyBlocks::corner";
    }

    return cor[i];
}

void ManyBlocks::set_corner(size_t i, Block* blk)
{
    if(i > 3) {
        throw "Illegal index for Block::cor in ManyBlocks::set_corner";
    }
     
     cor[i] = blk;
}

Block* ManyBlocks::back()
{
    return cor[1];
}

bool ManyBlocks::empty()
{
    return cor[0] == nullptr;
}

size_t ManyBlocks::size()
{
    return inv.size();
}

void ManyBlocks::link_block_corners()
{
    size_t N = inv.size();

    // forward sort by stop
    sort(2);
    for (size_t i = 0; i < N; i++) {
        inv[i]->cor[2] = (i == 0)     ? nullptr : inv[i - 1];
        inv[i]->cor[3] = (i == N - 1) ? nullptr : inv[i + 1];
    }
    // forward sort by start
    sort(0);
    for (size_t i = 0; i < N; i++) {
        inv[i]->cor[0] = (i == 0)     ? nullptr : inv[i - 1];
        inv[i]->cor[1] = (i == N - 1) ? nullptr : inv[i + 1];
    }
}

void ManyBlocks::link_corners(){
    size_t N = inv.size();
    try {
        for (size_t i = 0; i < 4; i++)
        {
            size_t k = i % 2 == 0 ? 0 : N - 1;
            Block* b = inv.at(k);
            while (b->corner(i) != nullptr)
            {
                b = b->corner(i);
            }
            cor[i] = b;
        }
    } catch (const std::out_of_range& e) {
        std::cerr << "Index error in " << __func__ << std::endl;
    }
}

void ManyBlocks::set_overlap_group(long &grpid)
{
    // Needed for determining overlaps and thus setids
    long maximum_stop = -1;
    // Loop through each Block in the linked list
    for (Block* blk = front(); blk != nullptr; blk = blk->next()) {
        // The stop position of the current interval
        long this_stop = blk->pos[1];
        // If the start is greater than the maximum stop, then the block is in
        // a new adjacency group. For this to work, Contig->block must be
        // sorted by start. This sort is performed in build_tree.
        if (blk->pos[0] > maximum_stop) {
            grpid++;
        }
        if (this_stop > maximum_stop) {
            maximum_stop = this_stop;
        }
        blk->grpid = grpid;
    }
}

    /** Link each node to its adjacent neighbors
     *
     * Link blocks to nearest non-overlapping up and downstream blocks
     *
     * For example, given these for blocks:
     *  |---a---|
     *            |--b--|
     *             |----c----|
     *                     |---d---|
     *                               |---e---|
     * a->adj := (nullptr, b)
     * b->adj := (a, e)
     * c->adj := (a, e)
     * d->adj := (a, e)
     * e->adj := (d, nullptr)
     */
void ManyBlocks::link_adjacent_blocks_directed(Direction d)
{
    // In diagrams:
    // <--- indicates a hi block
    // ---> indicates a lo block
    // All diagrams and comments relative to the d==HI direction

    if (cor[0] == nullptr || cor[1] == nullptr || cor[2] == nullptr || cor[3] == nullptr) {
        Rcpp::stop("Contig head must be set before link_adjacent_blocks is called\n");
    }

    // Transformed indices for Block->cor and Contig->cor
    int idx_a = (!d * 2) + !d; // - 0 previous/first element by start
    int idx_b = (!d * 2) +  d; // - 1 next/last element by start
    int idx_c = ( d * 2) + !d; // - 2 previous/first element by stop
    int idx_d = ( d * 2) +  d; // - 3 next/last element by stop

    Block *lo, *hi;

    lo = cor[idx_c]; // first element by stop
    hi = cor[idx_a]; // first element by start

    while (hi != nullptr) {

        //       --->
        // <---
        // OR
        // --->
        //   <---
        // This should should occur only at the beginning
        if (REL_LE(hi->pos[!d], lo->pos[d], d)) {
            hi->adj[!d] = nullptr;
            hi = hi->cor[idx_b];
        }

        //  lo     next
        // ---->  ---->
        //               <---
        // If next is closer, and not overlapping the hi, increment lo
        // You increment lo until it is adjacent to the current hi
        else if (REL_LT(lo->cor[idx_d]->pos[d], hi->pos[!d], d)) {
            lo = lo->cor[idx_d];
        }

        // --->
        //      <---
        // The current lo is next to, and not overlapping, current hi
        else {
            hi->adj[!d] = lo;
            hi = hi->cor[idx_b];
        }
    }
}

void ManyBlocks::link_adjacent_blocks()
{
    link_adjacent_blocks_directed(HI);
    link_adjacent_blocks_directed(LO);
}

void ManyBlocks::merge_overlaps()
{
    Block *lo, *hi;

    // iterate through all blocks
    for (lo = front(); lo != nullptr; lo = lo->next()) {
        // look ahead to find all doubly-overlapping blocks
        for (hi = lo->next(); hi != nullptr; hi = hi->next()) {
            if (!hi->overlap(lo)) {
                break;
            }
            if (hi->over->overlap(lo->over) && hi->over->parent == lo->over->parent) {
                Block::merge_block_a_into_b(hi, lo);
                hi = lo;
            }
        }
    }
}

void ManyBlocks::refresh()
{
    Block* first = inv[0];

    // NOTE: You cannot simply call front() to get the first defined block.

    // Get the first block that is defined
    for(auto &b : inv){
        if(b->over != nullptr){
            first = b;
            goto found;
        }
    }
throw "No legal blocks remain\n";
found:

    // Rewind (should be unnecessary, but I don't require the pool to be sorted)
    while(first->prev() != nullptr){
        first = first->prev();
    }

    // merging may invalidate the corners, so set to null
    cor = {{ nullptr }};

    // likewise, the Block array is invalidated, so clear this memory
    inv.clear();

    // refill it with the overlap-merged remaining blocks
    for(Block* b = first; b != nullptr; b = b->next()){
        // This works since the Blocks are zeroed in merge_block_a_into_b
        // And over must be defined in a well-formed block
        if(b->over != 0){
            inv.push_back(b);
        }
    }

    link_corners();
}
