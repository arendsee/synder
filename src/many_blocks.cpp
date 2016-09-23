#include "many_blocks.h"

ManyBlocks::ManyBlocks() { }

ManyBlocks::~ManyBlocks(){
    delete tree;
}

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

// TODO: set up a constant time alternative
size_t ManyBlocks::size()
{
    return inv.size();
}

void ManyBlocks::clear()
{
    inv.clear();
    delete tree;
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
    Block * b;
    size_t k;
    size_t N = inv.size();
    try {
        for (size_t i = 0; i < 4; i++)
        {
            k = i % 2 == 0 ? 0 : N - 1;
            b = inv.at(k);
            while (b->corner(i) != nullptr)
            {
                b = b->corner(i);
            }
            cor[i] = b;
        }
    } catch (const std::out_of_range& e) {
        cerr << "Index error in " << __func__ << endl;
    }
}

void ManyBlocks::set_overlap_group()
{
    // Holds current overlapping group id
    size_t grpid = 1;
    // Needed for determining overlaps and thus setids
    long maximum_stop = 0;
    // The stop position of the current interval
    long this_stop = 0;

    maximum_stop = 0;
    // Loop through each Block in the linked list
    for (Block* blk = front(); blk != nullptr; blk = blk->next()) {
        this_stop = blk->pos[1];
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
        fprintf(stderr, "Contig head must be set before link_adjacent_blocks is called\n");
        // fprintf(stderr, "genome=(%s) contig=(%s)\n", parent->parent->name.c_str(), parent->name.c_str());
        exit(EXIT_FAILURE);
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
            if (! hi->overlap(lo)) {
                break;
            }
            if (hi->over->overlap(lo->over) && hi->over->parent == lo->over->parent) {
                Block::merge_block_a_into_b(hi, lo);
                hi = lo;
            }
        }
    }

    size_t i = 0;
    for(; inv[i] == nullptr; i++){ }
    lo = inv[i];

    // rewind
    while(lo->prev() != nullptr){
        lo = lo->prev();
    }

    // clear pointer array
    inv.clear();

    // refill it with the overlap-merged remaining blocks
    while(lo != nullptr){
        inv.push_back(lo);    
        lo = lo->next();
    }

    // adjust corners as necessary
    link_corners();

}
