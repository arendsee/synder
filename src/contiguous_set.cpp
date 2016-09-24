#include "contiguous_set.h"

ContiguousSet::ContiguousSet() { }

/** Initialize a ContiguousSet from its first constituent
 *
 * Explicitly instantiating the ContiguousSet homolog is not really necessary,
 * since all information needed for its construction is present here.
 *
 */
ContiguousSet::ContiguousSet(Block* t_blk, size_t t_id)
    :
    LinkedInterval(
        t_blk->parent,
        0,
        t_blk->strand,
        t_id
    ),
    Interval(t_blk->pos[0], t_blk->pos[1]),
    size(1),
    ends({t_blk, t_blk})
{
    t_blk->cset = this;
}

ContiguousSet::ContiguousSet(ContiguousSet* cset)
    :
    LinkedInterval(
        cset->ends[0]->over->parent,
        cset->score,
        cset->ends[0]->over->strand,
        cset->id
    ),
    size(cset->size),
    ends({cset->ends[0]->over, cset->ends[1]->over})
{
    if(strand == '-'){
        std::swap(ends[0], ends[1]);
    }

    for(Block* blk = cset->ends[0]; blk != nullptr; blk = blk->cnr[1]){
        blk->over->cnr[0] = blk->cnr[0] == nullptr ? nullptr : blk->cnr[0]->over;
        blk->over->cnr[1] = blk->cnr[1] == nullptr ? nullptr : blk->cnr[1]->over;
        blk->over->cset = this;
    }

    pos[0] = ends[0]->pos[0];
    pos[1] = ends[1]->pos[1];
}

ContiguousSet::~ContiguousSet()
{
    if (prev != nullptr) {
        prev->next = next;
    }

    if (next != nullptr) {
        next->prev = prev;
    }
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
    return a->parent == b->parent &&
           (b->grpid - a->grpid) > (k+1);
}

bool ContiguousSet::are_contiguous(Block* blk_a, Block* blk_b, long k)
{
    long qdiff, tdiff, demerits;
    char ats, bts;
    long aqg, atg, bqg, btg;
    Feature *atc, *btc, *aqc, *bqc;

    // shortcuts to parents
    bqc = blk_b->parent;
    btc = blk_b->over->parent;
    aqc = blk_a->parent;
    atc = blk_a->over->parent;

    // shortcuts to strands
    bts = blk_b->over->strand;
    ats = blk_a->over->strand;

    // shortcuts to overlapping groups
    btg = blk_b->over->grpid;
    bqg = blk_b->grpid;
    atg = blk_a->over->grpid;
    aqg = blk_a->grpid;

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

    if (may_add) {
        force_add_block(blk_b);
    }
    return may_add;
}

void ContiguousSet::force_add_block(Block* blk_b)
{
    add_side_(blk_b);
    blk_b->cset = this;
    size++;
}

void ContiguousSet::add_side_(Block* b)
{
    assert(b->parent == parent);
    assert(b->strand == strand);

    Block* a;

    Direction d = static_cast<Direction> (b->strand == '+');

    if (b->pos[0] > pos[1]) {
        pos[1]     = b->pos[1];
        a          = ends[1];
        a->cnr[ d] = b;
        b->cnr[!d] = a;
        ends[1]    = b;
    } else if (b->pos[1] < pos[0]) {
        pos[0]     = b->pos[0];
        a          = ends[0];
        a->cnr[!d] = b;
        b->cnr[ d] = a;
        ends[0]    = b;
    } else {
        throw "Adding illegal or out-of-order block to ContiguousSet";
    }

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
    for (Block* x = a->corner(up); x != b; x = x->corner(up)) {
        if (x == nullptr) {
            throw "Foul magic in __func__:__LINE__";
        }
        if (
            ! x->overlap(a) &&
            ! x->overlap(b) &&
            x->over->parent == a->over->parent
        )
            return true;
    }
    return false;
}
