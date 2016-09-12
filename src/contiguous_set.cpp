#include "contiguous_set.h"

void add_block_to_ContiguousSet_side_(ContiguousSet * cset, Block * blk);

bool blocks_conflict(Block * a, Block * b);

ContiguousSet * init_ContiguousSet_side_(Block * blk);

void free_ContiguousSet(ContiguousSet * cset)
{
    if (cset->prev != NULL)
    {
        cset->prev->next = cset->next;
    }

    if (cset->next != NULL)
    {
        cset->next->prev = cset->prev;
    }

    if (cset->over != NULL)
    {
        cset->over->over = NULL;
        free_ContiguousSet(cset->over);
    }

    if (cset->parent->cset == cset)
    {
        cset->parent->cset = cset->prev != NULL ? cset->prev : cset->next;
    }

    cset->parent->clear_cset_tree();

    free(cset);
}

ContiguousSet * init_ContiguousSet_side_(Block * blk)
{
    if (blk != NULL)
    {
        ContiguousSet* cset = (ContiguousSet*)malloc(sizeof(ContiguousSet));
        cset->bounds[0]     = blk->pos[0];
        cset->bounds[1]     = blk->pos[1];
        cset->strand        = blk->strand;
        cset->parent        = blk->parent;
        cset->ends[0]       = blk;
        cset->ends[1]       = blk;
        cset->next          = NULL;
        cset->prev          = NULL;
        cset->size          = 1;
        cset->id            = 0;
        // add this ContiguousSet to the stack
        if (cset->parent->cset == NULL)
        {
            cset->parent->cset = cset;
        }
        else
        {
            cset->prev = cset->parent->cset;
            cset->parent->cset->next = cset;
            cset->parent->cset = cset;
        }
        blk->cset = cset;
        return cset;
    }
    else
    {
        fprintf(stderr, "WARNING: Cannot initialize ContiguousSet from NULL blk\n");
        return (ContiguousSet*) NULL;
    }
}
ContiguousSet * init_ContiguousSet(Block * blk)
{
    ContiguousSet * cset1 = init_ContiguousSet_side_(blk);
    ContiguousSet * cset2 = init_ContiguousSet_side_(blk->over);
    cset1->over           = cset2;
    cset2->over           = cset1;
    return cset1;
}

void print_ContiguousSet(ContiguousSet * cset)
{
    if (cset != NULL)
    {
        fprintf(
            stderr,
            "cset (n=%zu) - (%s, %zu, %zu, %c) <-> (%s, %zu, %zu, %c)\n",
            cset->size,
            cset->parent->name.c_str(),
            cset->bounds[0],
            cset->bounds[1],
            cset->ends[0]->strand,
            cset->over->parent->name.c_str(),
            cset->over->bounds[0],
            cset->over->bounds[1],
            cset->over->ends[0]->strand
        );
    }
}

bool strictly_forbidden(Block * a, Block * b, long k)
{
    return
        a->parent == b->parent &&
        ((long)b->grpid - (long)a->grpid) > (k+1);
}

bool are_contiguous(Block * blk_a, Block * blk_b, long k)
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

bool add_block_to_ContiguousSet(ContiguousSet * cset, Block * blk_b, long k)
{
    // the latermost element in the ContiguousSet
    Block * blk_a;
    // does blk_b meet the conditions for inclusion in the set?
    bool may_add;

    // It would be safer to search all elements in the set and find the one
    // closest to blk_b. But so long as we iterate through blocks that are
    // ordered, blk_b will always be nearest to cset->end[1].
    blk_a = cset->ends[1];

    // Determine if the blocks are contiguous. Eventually I may implement a few
    // distinct contiguity functions.
    may_add = are_contiguous(blk_a, blk_b, k);

    if (may_add)
    {
        // build query and target side sets
        add_block_to_ContiguousSet_side_(cset, blk_b);
        add_block_to_ContiguousSet_side_(cset->over, blk_b->over);
    }
    return may_add;
}



void add_block_to_ContiguousSet_side_(ContiguousSet * cset, Block * b)
{

    assert(b->parent == cset->parent);
    assert(b->strand == cset->strand);

    Block * a;

    Direction d = (Direction)(cset->strand == '+');

    if (b->pos[0] > cset->bounds[1])
    {
        cset->bounds[1] = b->pos[1];
        a               = cset->ends[1];
        a->cnr[ d]      = b;
        b->cnr[!d]      = a;
        cset->ends[1]   = b;
    }
    else if (b->pos[1] < cset->bounds[0])
    {
        cset->bounds[0] = b->pos[0];
        a               = cset->ends[0];
        a->cnr[!d]      = b;
        b->cnr[ d]      = a;
        cset->ends[0]   = b;
    }
    else
    {
        perror("Adding illegal or out-of-order block to ContiguousSet, dying\n");
        exit(EXIT_FAILURE);
    }

    b->cset = cset;
    cset->size++;
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
bool blocks_conflict(Block * a, Block * b)
{
    int up = (a->strand == '+') ? NEXT_START : PREV_STOP;
    Block * x = a->cor[up];
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
