#include "contiguous_set.h"

~ContiguousSet()
{
    if (prev != NULL)
    {
        prev->next = next;
    }

    if (next != NULL)
    {
        next->prev = prev;
    }

    if (over != NULL)
    {
        over->over = NULL;
        delete over;
    }

    if (parent->cset == this)
    {
        parent->cset = prev != NULL ? prev : next;
    }

    parent->clear_cset_tree();
}

ContiguousSet::ContiguousSet(Block* blk)
{
    if (blk != NULL)
    {
        pos[0]     = blk->pos[0];
        pos[1]     = blk->pos[1];
        strand     = blk->strand;
        parent     = blk->parent;
        ends[0]    = blk;
        ends[1]    = blk;
        next       = NULL;
        prev       = NULL;
        size       = 1;
        id         = 0;
        // add this ContiguousSet to the stack
        if (parent->cset == NULL)
        {
            parent->cset = this;
        }
        else
        {
            prev = parent->cset;
            parent->cset->next = this;
            parent->cset = this;
        }
        blk->cset = this;
    }
    else
    {
        fprintf(stderr, "WARNING: Cannot initialize ContiguousSet from NULL blk\n");
        exit(EXIT_FAILURE);
    }
}

ContiguousSet* ContiguousSet::make_pair(Block* blk)
{
    ContiguousSet* a = new ContiguousSet(blk);
    ContiguousSet* b = new ContiguousSet(blk->over);
    a->over = b;
    b->over = a;
    return a;
}

void ContiguousSet::print()
{
    fprintf(
        stderr,
        "cset (n=%zu) - (%s, %zu, %zu, %c) <-> (%s, %zu, %zu, %c)\n",
        size,
        parent->name.c_str(),
        pos[0],
        pos[1],
        ends[0]->strand,
        over->parent->name.c_str(),
        over->pos[0],
        over->pos[1],
        over->ends[0]->strand
    );
}

bool ContiguousSet::strictly_forbidden(Block* a, Block* b, long k)
{
    return
        a->parent == b->parent &&
        ((long)b->grpid - (long)a->grpid) > (k+1);
}

bool ContiguousSet::are_contiguous(Block* blk_a, Block* blk_b, long k)
{
    long qdiff, tdiff, demerits;
    char ats, bts;
    long aqg, atg, bqg, btg;
    Contig *atc, *btc, *aqc, *bqc;

    // interval variables for new query
    bqc =       blk_b->parent;
    bqg = (long)blk_b->grpid;
    // three interval variables for new target
    btc =       blk_b->over->parent;
    btg = (long)blk_b->over->grpid;
    bts =       blk_b->over->strand;

    // interval variables for previous query
    aqc =       blk_a->parent;
    aqg = (long)blk_a->grpid;
    // three interval variables for previously seen target
    atc =       blk_a->over->parent;
    atg = (long)blk_a->over->grpid;
    ats =       blk_a->over->strand;

    // qdiff and tdiff describe the adjacency of blocks relative to the
    // query are target contigs, respectively. Cases:
    // ---
    // diff <= -2 : blocks are not adjacent
    // diff == -1 : blocks are adjacent on reverse strand
    // diff ==  0 : blocks overlap
    // diff ==  1 : blocks are adjacent
    // diff >=  2 : blocks are not adjacent
    qdiff    = bqg - aqg;
    tdiff    = btg - atg;
    demerits = abs(tdiff) + qdiff - 2;

    return
        // non-overlapping
        qdiff    != 0   &&
        tdiff    != 0   &&
        // same target strand
        bts      == ats &&
        // same scaffolds
        aqc      == bqc &&
        atc      == btc &&
        // within an acceptable distance
        demerits <= k   &&
        (
            // going in the right direction
            (tdiff > 0 && bts == '+') ||
            (tdiff < 0 && bts == '-')
        ) &&
        // no cis jumpers
        ! blocks_conflict(blk_a->over, blk_b->over) &&
        ! blocks_conflict(blk_a, blk_b);
}

bool ContiguousSet::add_block(Block* blk_b, long k)
{
    // the latermost element in the ContiguousSet
    Block* blk_a;
    // does blk_b meet the conditions for inclusion in the set?
    bool may_add;

    // It would be safer to search all elements in the set and find the one
    // closest to blk_b. But so long as we iterate through blocks that are
    // ordered, blk_b will always be nearest to cset->end[1].
    blk_a = ends[1];

    // Determine if the blocks are contiguous. Eventually I may implement a few
    // distinct contiguity functions.
    may_add = are_contiguous(blk_a, blk_b, k);

    if (may_add)
    {
        // build query and target side sets
        add_block_to_ContiguousSet_side_(blk_b);
        add_block_to_ContiguousSet_side_(over, blk_b->over);
    }
    return may_add;
}



void ContiguousSet::add_block_side_(Block* b)
{

    assert(b->parent == parent);
    assert(b->strand == strand);

    Block* a;

    Direction d = (Direction)(strand == '+');

    if (b->pos[0] > pos[1])
    {
        pos[1]     = b->pos[1];
        a          = ends[1];
        a->cnr[ d] = b;
        b->cnr[!d] = a;
        ends[1]    = b;
    }
    else if (b->pos[1] < pos[0])
    {
        pos[0]     = b->pos[0];
        a          = ends[0];
        a->cnr[!d] = b;
        b->cnr[ d] = a;
        ends[0]    = b;
    }
    else
    {
        perror("Adding illegal or out-of-order block to ContiguousSet, dying\n");
        exit(EXIT_FAILURE);
    }

    b->cset = this;
    size++;
}

/* Determine if any non-overlapping target elements mapping to query exist
 * between TARGET blocks a and b
 *
 * If any interval maps to a region in )a,b(, and to any region on the query, return false
 *           a            z           b
 *     T  <=====>       <===>      <=====>
 *           |            |           |
 *           |            |           |
 *     Q  <=====>       <===>      <=====>
 *                  **Conflict!**
 *
 *           a            z           b
 *     T  <=====>       <===>      <=====>
 *           |            \___________|_________
 *           |                        |         \
 *     Q  <=====>                  <=====>     <===>
 *                                         **Conflict!**
 *
 *     Q2               <===>
 *           a            |           b
 *     T  <=====>       <===>      <=====>
 *           |            z           |
 *           |                        |
 *     Q  <=====>                  <=====>
 *                 **No Conflict**
 *
 * ???
 *           a            z           b
 *     T  <=====>       <===>      <=====>
 *           |             \_______   |
 *           |                     \  |
 *     Q  <=====>                 <==>--->
 *
 * ???
 *           a                        b
 *     T  <=====>                 <==>--->
 *           |             _______/   |
 *           |            /           |
 *     Q  <=====>      <===>      <======>
 *                       z
 *
 * ???
 *           a                        b
 *     T  <=====>                 <======>   <===>
 *           |             ___________|________/
 *           |            /           |
 *     Q  <=====>      <===>      <======>
 *                       z
 */
bool ContiguousSet::blocks_conflict(Block* a, Block* b)
{
    int up = (a->strand == '+') ? NEXT_START : PREV_STOP;
    Block* x = a->cor[up];
    for (; x != b; x = x->cor[up])
    {
        if (x == NULL)
        {
            fprintf(stderr, "Foul magic in __func__:__LINE__");
            exit(EXIT_FAILURE);
        }
        if (
            ! block_overlap(x, a) &&
            ! block_overlap(x, b) &&
            x->over->parent == a->over->parent
        )
            return true;
    }
    return false;
}
