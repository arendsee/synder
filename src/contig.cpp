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
ResultContig<T>* Contig::get_region(Bound& bound, IntervalTree<T>* tree, bool get_flanks)
{

    tree = build_tree(tree);

    // Search itree
    IntervalResult* res = tree->get_overlaps(&bound);

    // I want everything that overlaps the flanks if I am doing a Block search.
    // For ContiguousSet searches, I currently only want overlapping intervals.
    // TODO find a better solution
    if (get_flanks) {
        res = add_whatever_overlaps_flanks(res);
    }

    return res;
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



// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// TODO integrate all this into Contig
// TODO get rid of all the hand-rolled LLs

// TODO replace with Interval overlap_length (in SearchInterval class)
double overlap_length_ll{
    // If the intervals overlap
    if(a1 <= b2 && b1 <= a2) {
        // Find the lower bound of the overlapping region
        long a = a1 > b1 ? a1 : b1;
        // Find the upper bound of the overlapping region
        long b = a2 > b2 ? b2 : a2;
        // Return the overlapping interval length
        return b - a + 1;
    } else {
        return 0;
    }
}

// A local utility structure used filter and store contiguous sets
typedef struct CSList {
    struct CSList* next;
    Block* bound[2];
    ContiguousSet* cset;
} CSList;
CSList* init_empty_CSList();
CSList* init_CSList(Block* blk);
void add_blk_CSList(CSList* cslist, Block* blk);
void add_cset_CSList(CSList* cslist, ContiguousSet* cset, long bounds[2]);
void free_CSList(CSList* cslist);

typedef struct SI_Bound {
    long bound;
    int flag;
} SI_Bound;

SI_Bound* get_si_bound(
    long q,
    Block* blk_bounds[2],
    Direction d,
    bool inverted
);

int get_flag(SI_Bound* br[2]);

void find_search_intervals(Bound& bound, char* seqname)
{
    // Max and min blocks retrieved from itree.
    Block* blk_bounds[2];
    // Search interval boundary information
    SI_Bound* bound_results[2];
    // List of contiguous sets
    CSList *cslist;
    // A pointer to the root node of cslist (needed only for freeing the list)
    CSList *root;
    // Does starget strand == '-'?
    bool inverted;
    // Row output of itree
    ResultContig* rc;
    // Row output of ctree
    ResultContig* crc;
    // Search interval score
    double score;

    // Get blocks overlapping the query
    rc = qcon->get_region(bounds, qcon->block, true);

    // get list of highest and lowest members of each contiguous set
    cslist = init_empty_CSList();
    root = cslist;
    for(size_t i = 0; i < rc->size; i++) {
        add_blk_CSList(cslist, rc->block[i]);
    }

    // TODO what am I doing here? 
    crc = qcon->get_region(bounds, qcon->cset, false);
    if(! (crc->inbetween || crc->leftmost || crc->rightmost) ) {
        for(size_t i = 0; i < crc->size; i++) {
            add_cset_CSList(cslist, crc->cset[i], bounds);
        }
    }

    // Iterate through each contiguous set, for each find the search interval
    // For each contiguous set, there is exactly one search interval.
    for(; cslist != NULL; cslist = cslist->next) {

        blk_bounds[LO] = cslist->bound[LO];
        blk_bounds[HI] = cslist->bound[HI];

        inverted = blk_bounds[HI]->over->strand == '-';

        score = calculate_score(bounds, blk_bounds[LO]);

        bound_results[inverted ^ LO] =
            get_si_bound(bounds.pos[0], blk_bounds, LO, inverted);
        bound_results[inverted ^ HI] =
            get_si_bound(bounds.pos[1], blk_bounds, HI, inverted);


        printf("%s\t%s\t%zu\t%zu\t%s\t%zu\t%zu\t%c\t%lf\t%zu\t%i\t%i\t%i\n",
               // Output column ids:
               seqname,                                         //  1
               blk_bounds[LO]->parent->name.c_str(),            //  2
               bounds.pos[0] + Offsets::out_start,              //  3
               bounds.pos[1] + Offsets::out_stop,               //  4
               blk_bounds[LO]->over->parent->name.c_str(),      //  5
               bound_results[LO]->bound + Offsets::out_start,   //  6
               bound_results[HI]->bound + Offsets::out_stop,    //  7
               blk_bounds[LO]->over->strand,                    //  8
               score,                                           //  9
               blk_bounds[LO]->cset->id,                        // 10
               bound_results[LO]->flag,                         // 11
               bound_results[HI]->flag,                         // 12
               rc->inbetween || rc->leftmost || rc->rightmost   // 13
              );

        free(bound_results[0]);
        free(bound_results[1]);

    }

    free_ResultContig(rc);
    free_ResultContig(crc);
    free_CSList(root);
}

SI_Bound* init_SI_Bound(long bound, int flag)
{
    SI_Bound* br = (SI_Bound *)malloc(sizeof(SI_Bound));
    br->bound = bound;
    br->flag = flag;
    return br;
}

SI_Bound* get_si_bound(
    long q,
    Block* blk_bounds[2],
    Direction d,
    bool inverted)
{
    // Invert orientation mapping to target if search interval is inverted
    Direction vd = (Direction) (inverted ? !d : d);
    // See contiguous.h
    int flag = 0;
    // non-zero to ease debugging
    long bound = 444444;

    long set_bounds[2];
    set_bounds[0] = blk_bounds[0]->cset->bounds[0];
    set_bounds[1] = blk_bounds[0]->cset->bounds[1];

    // All diagrams are shown for the d=HI case, take the mirror image fr d=LO.
    //
    // KEY:
    // |x--  --y| - bounds of the contiguous block; start==x, stop==y
    // a========b - a syntenic block with start == a and stop == b
    //   <---q    - the query interval, with stop == q (start doesn't matter)
    // a==b--c==d - query bounding blocks in the contiguous set
    // [===]      - a non-bounding block in the same contiguous set
    //  ...  F=== - nearest non-adjacent block ***ON THE TARGET***, F=start
    //      ^     - search interval bound
    //
    // Possible snap positions (relative to query)
    //   |x...[===]-----a=======b-----c=======d-----[===]...y|  ...  F===
    //                 ^        ^    ^        ^    ^                ^

    // Positions of a, b, c, and d (as shown above)
    long pnt_a = blk_bounds[!d]->pos[!d];
    long pnt_b = blk_bounds[!d]->pos[ d];
    long pnt_c = blk_bounds[ d]->pos[!d];
    long pnt_d = blk_bounds[ d]->pos[ d];


    // This may occur when there is only one element in the ContiguousSet
    //          |x-----y|
    //       ...a=======b
    //         ^
    //   <---q
    // q < x
    if(REL_LT(q, set_bounds[!d], d)) {
        bound = blk_bounds[!d]->over->pos[!vd];
        flag = UNBOUND;
    }

    //   |x---[===]-------a=======b...y|   ...    F===
    //                   ^
    //              --q
    // q < a
    else if(REL_LT(q, pnt_a, d)) {
        bound = blk_bounds[!d]->over->pos[!vd];
        flag = BOUND;
    }

    //   |x...a=======b...y|   ...    F===
    //                ^
    //        <---q
    // q < b
    else if(REL_LE(q, pnt_b, d)) {
        bound = blk_bounds[!d]->over->pos[vd];
        flag = ANCHORED;
    }

    //   |x...a=======b-------c=======d...y|  ...  F===
    //                       ^
    //                  <--q
    // q < c && q > b
    //   (q > b test required since blk_bounds[LO] can equal blk_bounds[HI])
    else if(REL_LT(q, pnt_c, d) && REL_GT(q, pnt_b, d)) {
        bound = blk_bounds[d]->over->pos[!vd];
        flag = BOUND;
    }

    //   |x...a=======b-------c=======d...y|  ...  F===
    //                                ^
    //             <--------------q
    // q < d
    else if(REL_LE(q, pnt_d, d)) {
        bound = blk_bounds[d]->over->pos[vd];
        flag = ANCHORED;
    }

    //   |x...a=======b-------c=======d-------[===]...y|  ...  F===
    //                                       ^
    //              <----------------------q
    // q < y, (which implies there is a node after d)
    else if(REL_LE(q, set_bounds[d], d)) {
        bound = blk_bounds[d]->cnr[d]->over->pos[!vd];
        flag = BOUND;
    }

    // If none of the above, the bound beyond anything in the contiguous set
    // In this case, the hi and lo Contiguous nodes will be the same
    else {
        // Get nearest non-overlapping sequence
        Block* downstream_blk = blk_bounds[d]->over->adj[vd];

        // adjacent block on TARGET side exists
        //    |x...--a=======b|
        //    |x...--c=======d|  ...  F===
        //                           ^
        //                    <---q
        if(downstream_blk != NULL) {
            flag = UNBOUND;
            bound = downstream_blk->pos[!vd];
        }
        //    |x...--a=======b|
        //    |x...--c=======d|  ...  THE_END
        //                    <---q
        // query is further out than ANYTHING in the synteny map
        else {
            bound = vd ? blk_bounds[d]->over->parent->length - 1 : 0;
            flag = blk_bounds[d]->over->pos[vd] == bound ? BEYOND : EXTREME;
        }
    }

    return init_SI_Bound(bound, flag);
}


// Arguments:
//  near - the distance from the block to the expected near end of the query
//  far  - the distance from the block to the expected far end of the query
//               b1        b2  |  b1        b2
//      a1  a2   |=========|   |  |=========|    a1  a2
//      |===|    |             |                 |===|
//      |<------>| far         |       near |<-->|
//          |<-->| near        |        far |<------>|
double _flank_area(long near, long far, double k)
{

    // If far <= 0, this means there is no interval to score in this direction
    if(far > 0) {
        // Adjust the near boundary, if needed, for example:
        //             b1        b2
        //             |=========|
        //      a1       |   a2
        //      |============|
        //      |<------>|     far
        //               |<->| near (negative value)
        //               |     snap near to relative 0
        near = near > 0 ? near : 0;

        // the weight falls exponentially with distance from the query, e.g.
        // $$ \int_{near}^{far} exp(-kx) dx $$
        // which evaluates to the following:
        return (1 / k) * ( exp(-1 * k * near) - exp(-1 * k * far) );
    }
    return 0;
}

double calculate_score(Bound& bound, Block* blk)
{

    long a1 = bound.pos[0];
    long a2 = bound.pos[1];

    double weighted_length;
    long b1, b2, actual_length;
    double k = 0.001;
    double score = 0;

    if(blk == NULL)
        return score;

    // rewind
    while(blk->cnr[0] != NULL) {
        blk = blk->cnr[0];
    }

    for(; blk != NULL ; blk = blk->cnr[1]) {
        b1 = blk->pos[0];
        b2 = blk->pos[1];

        //               a1        a2
        //      b1  b2   |=========|    query interval
        //      |===|    |              syntenic interval
        //      |---|....|              interval to score
        // near difference := i1 = a1 - b2
        // far difference  := i2 = a1 - b1
        // i1 and i2 may be negative, if query is not in the above position
        weighted_length = _flank_area(a1 - b2, a1 - b1, k) +

                          //    a1        a2
                          //    |=========|    b1  b2     query interval
                          //              |....|---|      syntenic interval
                          //                   |===|      interval to score
                          // near difference := i1 = b1 - a2
                          // far difference  := i2 = b1 - a1
                          _flank_area(b1 - a1, b2 - a1, k) +

                          //         a1        a2
                          //         |=========|         query interval
                          //             |=====|=====|   syntenic interval
                          //             b1    |     b1
                          //             |-----|         overlapping interval
                          //             i1    i2
                          overlap_length_ll(a1, a2, b1, b2);

        // NOTE: I am kind of adding length to area here, but it actually
        // works. `_flank_area` returns the area of a segment of the base
        // exponentional (i.e. where f(0) = 1). Multiplying this base
        // exponential area by the syntenic link score, gives the final score
        // for the non-overlapping segment. In the same way, multiplying the
        // overlapping segment length by the score, gives the overlapping
        // segmental score.

        actual_length = blk->pos[1] - blk->pos[0] + 1;

        score += blk->score * weighted_length / actual_length;
    }
    return score;
}


CSList* init_empty_CSList()
{
    CSList* cslist   = (CSList *)malloc(sizeof(CSList));
    cslist->next      = NULL;
    cslist->bound[LO] = NULL;
    cslist->bound[HI] = NULL;
    cslist->cset      = NULL;
    return(cslist);
}

CSList* init_CSList(Block* blk)
{
    CSList* cslist = (CSList *)malloc(sizeof(CSList));
    cslist->next = NULL;
    cslist->bound[LO] = blk;
    cslist->bound[HI] = blk;
    cslist->cset = blk->cset;
    return(cslist);
}

void add_cset_CSList(CSList* cslist, ContiguousSet* cset, long bounds[2])
{
    if(cset == NULL) {
        // do nothing
    } else if(cset == cslist->cset) {
        // do nothing
    } else if(cslist->next == NULL) {
        cslist->next = init_empty_CSList();
        cslist       = cslist->next;
        cslist->cset = cset;

        long lo = LONG_MIN;
        long hi = LONG_MAX;
        Block* blk = cset->ends[0];
        // blk should never be NULL
        assert(blk != NULL);
        // if ContiguousSet is valid, ends[0] should have no prev
        assert(blk->cnr[0] == NULL);
        for(; blk != NULL; blk = blk->cnr[1]) {
            if(blk->pos[0] <= bounds[0] && blk->pos[0] >= lo) {
                cslist->bound[LO] = blk;
                lo = blk->pos[0];
            }
            if(blk->pos[1] >= bounds[1] && blk->pos[1] <= hi) {
                cslist->bound[HI] = blk;
                hi = blk->pos[1];
            }
        }
    } else {
        add_cset_CSList(cslist->next, cset, bounds);
    }
}

void add_blk_CSList(CSList* cslist, Block* blk)
{
    // first entry of empty CSList
    if(blk == NULL) {
        fprintf(stderr, "I got a null block\n");
        exit(EXIT_FAILURE);
    } else if(cslist->bound[LO] == NULL && cslist->bound[HI] == NULL) {
        cslist->bound[LO] = blk;
        cslist->bound[HI] = blk;
        cslist->cset = blk->cset;
    } else if(cslist->cset == blk->cset) {
        if(cslist->bound[HI] == NULL || blk->pos[LO] > cslist->bound[HI]->pos[HI]) {
            cslist->bound[HI] = blk;
        }
        if(cslist->bound[LO] == NULL || blk->pos[HI] < cslist->bound[LO]->pos[LO]) {
            cslist->bound[LO] = blk;
        }
    } else if(cslist->next == NULL) {
        cslist->next = init_CSList(blk);
    } else {
        add_blk_CSList(cslist->next, blk);
    }
}

void free_CSList(CSList* cslist)
{
    if(cslist->next != NULL)
        free_CSList(cslist->next);
    free(cslist);
}
