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

Contig::Contig(std::string new_name, Genome* new_parent)
    : name(new_name), parent(new_parent)
{
    length = default_length;
    itree  = NULL;
    ctree  = NULL;
    cset   = NULL;
}

Contig::~Contig()
{
    while (cset != NULL) {
        free_ContiguousSet(cset);
    }
    if (itree != NULL) {
        delete itree;
    }
    if (ctree != NULL) {
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
    for (; c != NULL; c = c->next) {
        fprintf(stderr, "  -- ");
        print_ContiguousSet(c);
    }

    if (print_blocks) {
        int d = forward ? 1 : 3; // next by start or next by stop
        Block* blk = cor[d - 1]; // prev by start or prev by stop
        for (; blk != NULL; blk = blk->cor[d]) {
            print_Block(blk);
        }
    }
}

template<typename T>
T* add_whatever_overlaps_flanks(T* res)
{
    // itree returns the flanks for queries that overlap nothing. However, I
    // need all the intervals that overlap these flanks as well.
    T* tmp_a = NULL;
    T* tmp_b = NULL;

    if (res->inbetween) {
        // If inbetween, itree should have returned the two flanking blocks
        if (res->iv.size() == 2) {
            tmp_a = res->tree->get_overlaps(res->iv[0]);
            tmp_b = res->tree->get_overlaps(res->iv[1]);
            res->iv.clear();
            res->iv = tmp_a->iv;
            res->iv.insert(
                res->iv.end(),
                tmp_b->iv.begin(),
                tmp_b->iv.end()
            );
        } else {
            fprintf(stderr, "itree is broken, should return exactly 2 intervals for inbetween cases\n");
            exit(EXIT_FAILURE);
        }
    } else if (res->leftmost || res->rightmost) {
        if (res->iv.size() == 1) {
            tmp_a = res->tree->get_overlaps(res->iv[0]);
            res->iv = tmp_a->iv;
        } else {
            fprintf(stderr, "itree is broken, should return only 1 interval for left/rightmost cases\n");
            exit(EXIT_FAILURE);
        }
    }

    if (tmp_a != NULL)
        delete tmp_a;

    if (tmp_b != NULL)
        delete tmp_b;

    return res;
}

template<typename T>
void Contig::build_block_itree()
{
    if (itree == NULL) {

        std::vector<Block*> intervals;

        for (Block* blk = cor[0]; blk != NULL; blk = blk->cor[1]) {
            intervals.push_back(blk);
        }

        itree = new IntervalTree(intervals);
    }
}

void Contig::build_itree()
{
    if (ctree == NULL) {

        std::vector<ContiguousSet*> intervals;

        for (ContiguousSet* c = cset; c != NULL; c = c->next) {
            intervals.push_back(c);
        }

        ctree = new IntervalTree(intervals);
    }
}

template <class T>
ResultContig<T>* Contig::get_region(long a, long b, IntervalTree<T>* tree)
{

    tree = build_tree(tree);

    // Search itree
    Bound inv(a, b);
    IntervalResult* res = tree->get_overlaps(&inv);

    // I want everything that overlaps the flanks if I am doing a Block search.
    // For ContiguousSet searches, I currently only want overlapping intervals.
    // TODO find a better solution
    if (!is_cset) {
        res = add_whatever_overlaps_flanks(res);
    }

    return res;
}

IntervalResult<ContiguousSet>* get_cset_region(long a, long b){
    get_region(a, b);
}

IntervalResult<Block>* get_block_region(long a, long b){
    get_region(a, b);
}

void Contig::link_contiguous_blocks(long k, size_t &setid)
{

    std::list<ContiguousSet*> csets;
    std::list<ContiguousSet*>::iterator iter = csets.begin();

    for (Block* blk = cor[0]; blk != NULL; blk = blk->cor[1]) {
        iter = csets.begin();
        while (true) {
            if (iter == csets.end()) {
                // if block fits in no set, create a new one
                csets.push_front(new ContiguousSet(blk));
                break;
            }
            // if block has joined a set
            else if ((*iter)->add_block(blk, k)) {
                break;
            }
            // if set terminates
            else if (strictly_forbidden((*iter)->ends[1], blk, k)) {
                iter = csets.erase(iter);
            } else {
                iter++;
            }
        }
    }

    ContiguousSet* cset_ptr = *csets.begin();
    // rewind - TODO - build the csets such that this isn't necessary
    while (cset_ptr->prev != NULL) {
        cset_ptr = cset_ptr->prev;
    }
    cset = cset_ptr;
    while (cset_ptr != NULL) {
        setid++;
        cset_ptr->id = setid;
        cset_ptr->over->id = setid;
        cset_ptr = cset_ptr->next;
    }
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
    for (lo = cor[0]; lo != NULL; lo = lo->cor[1]) {
        // look ahead to find all doubly-overlapping blocks
        for (hi = lo->cor[1]; hi != NULL; hi = hi->cor[1]) {
            if (! hi->overlap(lo)) {
                break;
            }
            if (hi->over->overlap(lo->over) && hi->over->parent == lo->over->parent) {
                merge_block_a_into_b(hi, lo);
                hi = lo;
            }
        }
    }
}

void Contig::sort_blocks(bool by_stop)
{
    auto cmp = by_stop ? Block::cmp_start : Block::cmp_stop;
    std::sort(blocks.begin(), blocks.end(), cmp);
}

void Contig::clear_cset_tree()
{
    if (ctree != NULL) {
        delete ctree;
        ctree = NULL;
    }
}
