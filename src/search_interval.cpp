#include "search_interval.h"

void SearchInterval::SearchInterval(const Feature& t_feat)
    : feat(t_feat)
{ }

void build_search_interval(Bound& bound, Block* blk_bound[2], char* seqname, bool new_inbetween)
{
    // Search interval boundary information
    SI_Bound* bound_results[2];

    bool inverted = blk_bounds[HI]->over->strand == '-';

    float score = calculate_score(bounds, blk_bounds[LO]);

    bound_results[inverted ^ LO] =
        get_si_bound(bounds.pos[LO], blk_bounds, LO, inverted);
    bound_results[inverted ^ HI] =
        get_si_bound(bounds.pos[HI], blk_bounds, HI, inverted);

    qbound = &bound;
    cset = blk_bounds[LO]->cset;
    inbetween = new_inbetween;
    pos[0]  = bound_results[0]->bound;
    pos[1]  = bound_results[1]->bound;
    flag[0] = bound_result[0]->flag;
    flag[1] = bound_result[1]->flag;
}

void SearchInterval::print()
{
    printf("%s\t%s\t%zu\t%zu\t%s\t%zu\t%zu\t%c\t%lf\t%zu\t%i\t%i\t%i\n",
           // Output column ids:
           feat.feature_name,                             //  1
           feat.contig_name,          //  2
           qbound->pos[0] + Offsets::out_start, //  3
           qbound->pos[1] + Offsets::out_stop,  //  4
           cset->over->parent->name.c_str(),    //  5
           pos[0] + Offsets::out_start,         //  6
           pos[1] + Offsets::out_stop,          //  7
           cset->over->strand,                  //  8
           score,                               //  9
           cset->id,                            // 10
           flag[0],                             // 11
           flag[1],                             // 12
           inbetween                            // 13 TODO fix the smell
          );
}

SI_Bound* SearchInterval::init_SI_Bound(long bound, int flag)
{
    SI_Bound* br = (SI_Bound *)malloc(sizeof(SI_Bound));
    br->bound = bound;
    br->flag = flag;
    return br;
}

SI_Bound* SearchInterval::get_si_bound(
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
double SearchInterval::flank_area(long near, long far, double k)
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

double SearchInterval::calculate_score(Bound& bound, Block* blk)
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
        weighted_length = flank_area(a1 - b2, a1 - b1, k) +

                          //    a1        a2
                          //    |=========|    b1  b2     query interval
                          //              |....|---|      syntenic interval
                          //                   |===|      interval to score
                          // near difference := i1 = b1 - a2
                          // far difference  := i2 = b1 - a1
                          flank_area(b1 - a1, b2 - a1, k) +

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

CSList* SearchInterval::init_empty_CSList()
{
    CSList* cslist   = (CSList *)malloc(sizeof(CSList));
    cslist->next      = NULL;
    cslist->bound[LO] = NULL;
    cslist->bound[HI] = NULL;
    cslist->cset      = NULL;
    return(cslist);
}

CSList* SearchInterval::init_CSList(Block* blk)
{
    CSList* cslist = (CSList *)malloc(sizeof(CSList));
    cslist->next = NULL;
    cslist->bound[LO] = blk;
    cslist->bound[HI] = blk;
    cslist->cset = blk->cset;
    return(cslist);
}

void SearchInterval::add_cset_CSList(CSList* cslist, ContiguousSet* cset, long bounds[2])
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

void SearchInterval::add_blk_CSList(CSList* cslist, Block* blk)
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

void SearchInterval::free_CSList(CSList* cslist)
{
    if(cslist->next != NULL)
        free_CSList(cslist->next);
    free(cslist);
}
