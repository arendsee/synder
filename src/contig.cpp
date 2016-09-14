#include "contig.h"

Contig::Contig()
{
    length = default_length;
    name   = "";
    parent = NULL;
    itree  = NULL;
    ctree  = NULL;
    cset   = NULL;
}

Contig::Contig(std::string new_name, Genome * new_parent)
    : name(new_name), parent(new_parent)
{
    length = default_length;
    itree  = NULL;
    ctree  = NULL;
    cset   = NULL;
}

Contig::~Contig()
{
    while (cset != NULL)
    {
        free_ContiguousSet(cset);
    }
    if (itree != NULL)
    {
        delete itree;
    }
    if (ctree != NULL)
    {
        delete ctree;
    }
}

void Contig::print(bool forward, bool print_blocks)
{
    fprintf(
        stderr,
        "$ %s size=%lu length=%lu cor=[%zu,%zu,%zu,%zu]\n",
        name.c_str(),
        block.size(),
        length,
        cor[0]->linkid,
        cor[1]->linkid,
        cor[2]->linkid,
        cor[3]->linkid
    );

    ContiguousSet* c = cset;
    for (; c != NULL; c = c->next)
    {
        fprintf(stderr, "  -- ");
        print_ContiguousSet(c);
    }

    if (print_blocks)
    {
        int d = forward ? 1 : 3; // next by start or next by stop
        Block* blk = cor[d - 1]; // prev by start or prev by stop
        for (; blk != NULL; blk = blk->cor[d])
        {
            print_Block(blk);
        }
    }
}

void Contig::build_block_itree()
{
    // build itree if necessary
    if (itree == NULL)
    {
        Block* blk = cor[0];
        // count the number of blks in the linked list
        // since blocks can be deleted, we cannot just use con->size
        size_t n = 0;
        for (; blk != NULL; blk = blk->cor[1])
        {
            n++;
        }

        // allocate Interval memory pool for IntervalTree
        // IntervalTree's destructor will free this block
        Interval* intervals = (Interval*)malloc(n * sizeof(Interval));

        // index for array of intervals
        size_t i = 0;
        // reset blk, since was wound to NULL above
        blk = cor[0];
        // map the Block list to the new Interval pool
        for (; blk != NULL; blk = blk->cor[1], i++)
        {
            intervals[i].start = blk->pos[0];
            intervals[i].stop  = blk->pos[1];
            intervals[i].link  = (void*)blk;
        }

        itree = new IntervalTree(intervals, n);
    }
}

void Contig::build_cset_itree()
{
    // Build ctree if necessary
    if (ctree == NULL)
    {
        ContiguousSet * c = cset;
        // rewind cset if necessary
        while (c->prev != NULL)
        {
            c = c->prev;
        }
        size_t n = 0;
        for (; c != NULL; c = c->next)
        {
            n++;
        }
        // Reset cset, since was wound to NULL above
        c = cset;

        // Map the Block list to the new IA structure
        Interval* intervals = (Interval*)malloc(n * sizeof(Interval));

        for (size_t i = 0; c != NULL; c = c->next, i++)
        {
            intervals[i].start = c->bounds[0];
            intervals[i].stop  = c->bounds[1];
            intervals[i].link  = (void*)c;
        }

        ctree = new IntervalTree(intervals, n);
    }
}

void Contig::link_contiguous_blocks(long k, size_t &setid)
{

    std::list<ContiguousSet*> csets;
    std::list<ContiguousSet*>::iterator iter = csets.begin();

    for (Block * blk = cor[0]; blk != NULL; blk = blk->cor[1])
    {
        iter = csets.begin();
        while (true)
        {
            if (iter == csets.end())
            {
                // if block fits in no set, create a new one
                csets.push_front(init_ContiguousSet(blk));
                break;
            }
            // if block has joined a set
            else if (add_block_to_ContiguousSet(*iter, blk, k))
            {
                break;
            }
            // if set terminates
            else if (strictly_forbidden((*iter)->ends[1], blk, k))
            {
                iter = csets.erase(iter);
            }
            else
            {
                iter++;
            }
        }
    }

    ContiguousSet* cset_ptr = *csets.begin();
    // rewind - TODO - build the csets so this isn't necessary
    while (cset_ptr->prev != NULL)
    {
        cset_ptr = cset_ptr->prev;
    }
    cset = cset_ptr;
    while (cset_ptr != NULL)
    {
        setid++;
        cset_ptr->id = setid;
        cset_ptr->over->id = setid;
        cset_ptr = cset_ptr->next;
    }
}

ResultContig* Contig::get_region(long a, long b, bool is_cset)
{

    IntervalTree * tree;

    if (is_cset)
    {
        build_cset_itree();
        tree = ctree;
    }
    else
    {
        build_block_itree();
        tree = itree;
    }

    // Search itree
    Interval inv(a, b);
    IntervalResult* res = tree->get_overlaps(&inv);

    // itree returns the flanks for queries that overlap nothing. However, I
    // need all the intervals that overlap these flanks as well.
    IntervalResult* tmp_a = NULL;
    IntervalResult* tmp_b = NULL;

    // I want everything that overlaps the flanks if I am doing a Block search.
    // For ContiguousSet searches, I currently only want overlapping intervals.
    if (!is_cset)
    {
        if (res->inbetween)
        {
            // If inbetween, itree should have returned the two flanking blocks
            if (res->iv.size() == 2)
            {
                tmp_a = itree->get_overlaps(res->iv[0]);
                tmp_b = itree->get_overlaps(res->iv[1]);
                res->iv.clear();
                res->iv = tmp_a->iv;
                res->iv.insert(
                    res->iv.end(),
                    tmp_b->iv.begin(),
                    tmp_b->iv.end()
                );
            }
            else
            {
                fprintf(stderr, "itree is broken, should return exactly 2 intervals for inbetween cases\n");
                exit(EXIT_FAILURE);
            }
        }
        else if (res->leftmost || res->rightmost)
        {
            if (res->iv.size() == 1)
            {
                tmp_a = itree->get_overlaps(res->iv[0]);
                res->iv = tmp_a->iv;
            }
            else
            {
                fprintf(stderr, "itree is broken, should return only 1 interval for left/rightmost cases\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    ResultContig* resultcontig = init_ResultContig(this, res, is_cset);

    if (tmp_a != NULL)
        delete tmp_a;
    if (tmp_b != NULL)
        delete tmp_b;

    delete res;

    return resultcontig;
}

long Contig::count_overlaps(long a, long b)
{
    build_block_itree();
    Interval inv(a, b);
    long count = itree->count_overlaps(&inv);
    return count;
}


void Contig::merge_doubly_overlapping_blocks()
{
    Block *lo, *hi;

    // iterate through all blocks
    for (lo = cor[0]; lo != NULL; lo = lo->cor[1])
    {
        // look ahead to find all doubly-overlapping blocks
        for (hi = lo->cor[1]; hi != NULL; hi = hi->cor[1])
        {
            if (! block_overlap(hi, lo))
            {
                break;
            }
            if (block_overlap(hi->over, lo->over) && hi->over->parent == lo->over->parent)
            {
                merge_block_a_into_b(hi, lo);
                hi = lo;
            }
        }
    }
}

void Contig::sort_blocks(Block** blocks, size_t size, bool by_stop)
{
    if (blocks != NULL)
    {
        if (by_stop)
        {
            qsort(blocks, size, sizeof(Block*), block_cmp_stop);
        }
        else
        {
            qsort(blocks, size, sizeof(Block*), block_cmp_start);
        }
    }
}

void Contig::clear_cset_tree()
{
    if (ctree != NULL)
    {
        delete ctree;
        ctree = NULL;
    }
}
