#include "contig.h"

Contig::Contig()
{
    length = default_length;
    name   = "";
    parent = NULL;
    // TODO handle IntervalSet initialization
}

Contig::Contig(std::string new_name, Genome* new_parent)
    : parent(new_parent), name(new_name)
{
    length = default_length;
    // TODO handle IntervalSet initialization
}

Contig::~Contig()
{
    // TODO I don't think there actually is anything to remove
}

void Contig::print(bool forward, bool print_blocks)
{
    // fprintf(
    //     stderr,
    //     "$ %s size=%lu length=%lu cor=[%zu,%zu,%zu,%zu]\n",
    //     name.c_str(),
    //     block.size(),
    //     length,
    //     cor[0]->linkid,
    //     cor[1]->linkid,
    //     cor[2]->linkid,
    //     cor[3]->linkid
    // );
    // 
    // ContiguousSet* c = cset;
    // for (; c != NULL; c = c->next) {
    //     fprintf(stderr, "  -- ");
    //     print_ContiguousSet(c);
    // }
    // 
    // if (print_blocks) {
    //     int d = forward ? 1 : 3; // next by start or next by stop
    //     Block* blk = cor[d - 1]; // prev by start or prev by stop
    //     for (; blk != NULL; blk = blk->cor[d]) {
    //         print_Block(blk);
    //     }
    // }
}

void Contig::link_contiguous_blocks(long k, size_t &setid)
{
    // 
    // std::list<ContiguousSet*> csets;
    // std::list<ContiguousSet*>::iterator iter = csets.begin();
    // 
    // for (Block* blk = cor[0]; blk != NULL; blk = blk->cor[1]) {
    //     iter = csets.begin();
    //     while (true) {
    //         if (iter == csets.end()) {
    //             // if block fits in no set, create a new one
    //             csets.push_front(new ContiguousSet(blk));
    //             break;
    //         }
    //         // if block has joined a set
    //         else if ((*iter)->add_block(blk, k)) {
    //             break;
    //         }
    //         // if set terminates
    //         else if (strictly_forbidden((*iter)->ends[1], blk, k)) {
    //             iter = csets.erase(iter);
    //         } else {
    //             iter++;
    //         }
    //     }
    // }
    // 
    // ContiguousSet* cset_ptr = *csets.begin();
    // // rewind - TODO - build the csets such that this isn't necessary
    // while (cset_ptr->prev != NULL) {
    //     cset_ptr = cset_ptr->prev;
    // }
    // cset = cset_ptr;
    // while (cset_ptr != NULL) {
    //     setid++;
    //     cset_ptr->id = setid;
    //     cset_ptr->over->id = setid;
    //     cset_ptr = cset_ptr->next;
    // }
}

void Contig::merge_doubly_overlapping_blocks()
{
    // Block *lo, *hi;
    // 
    // // iterate through all blocks
    // for (lo = cor[0]; lo != NULL; lo = lo->cor[1]) {
    //     // look ahead to find all doubly-overlapping blocks
    //     for (hi = lo->cor[1]; hi != NULL; hi = hi->cor[1]) {
    //         if (! hi->overlap(lo)) {
    //             break;
    //         }
    //         if (hi->over->overlap(lo->over) && hi->over->parent == lo->over->parent) {
    //             merge_block_a_into_b(hi, lo);
    //             hi = lo;
    //         }
    //     }
    // }
}

void Contig::sort_blocks(bool by_stop)
{
    // auto cmp = by_stop ? Block::cmp_start : Block::cmp_stop;
    // std::sort(blocks.begin(), blocks.end(), cmp);
}

void Contig::clear_cset_tree()
{
    // TODO implement
}

void Contig::count(Bound& bound, char* seqname)
{
    // build_block_itree();
    // long count = itree->count_overlaps(bound);
    // printf("%s\t%zu\n", seqname, count);
}

void Contig::map(Bound& bound, char* seqname)
{
    // Block *qblk, *tblk;
    // ContigResult* rc = get_region(bound, false);
    // bool missing = rc->inbetween || rc->leftmost || rc->rightmost;
    // 
    // for (size_t i = 0; i < rc->size; i++) {
    //     qblk = rc->block[i];
    //     if (qblk != NULL) {
    //         tblk = qblk->over;
    //         printf("%s %s %zu %zu %d\n",
    //                seqname,
    //                tblk->parent->name.c_str(),
    //                tblk->pos[0] + Offsets::out_start,
    //                tblk->pos[1] + Offsets::out_stop,
    //                missing
    //               );
    //     }
    // }
    // free(rc);
}

void Contig::find_search_intervals(Bound& bounds, char* seqname){
    // // List of contiguous sets
    // CSList *cslist;
    // // A pointer to the root node of cslist (needed only for freeing the list)
    // CSList *root;
    // // Does starget strand == '-'?
    // bool inverted;
    // // Row output of itree
    // ResultContig* rc;
    // // Row output of ctree
    // ResultContig* crc;
    // // Output search interval
    // SearchInterval si(this);
    // 
    // // Get blocks overlapping the query
    // rc = get_region(bounds, block, true);
    // 
    // // get list of highest and lowest members of each contiguous set
    // cslist = init_empty_CSList();
    // root = cslist;
    // for(size_t i = 0; i < rc->size; i++) {
    //     add_blk_CSList(cslist, rc->block[i]);
    // }
    // 
    // // TODO what am I doing here?
    // // Merge all this crap into the SearchInterval class
    // crc = get_region(bounds, cset, false);
    // if(! (crc->inbetween || crc->leftmost || crc->rightmost) ) {
    //     for(size_t i = 0; i < crc->size; i++) {
    //         add_cset_CSList(cslist, crc->cset[i], bounds);
    //     }
    // }
    // 
    // // Iterate through each contiguous set, for each find the search interval
    // // For each contiguous set, there is exactly one search interval, or into a new SearchIntervalSet class
    // for(; cslist != NULL; cslist = cslist->next) {
    // 
    //     blk_bounds[LO] = cslist->bound[LO];
    //     blk_bounds[HI] = cslist->bound[HI];
    // 
    //     // TODO fix the smell
    //     si = build_search_interval(bounds, blk_bounds, seqname, rc->inbetween || rc->leftmost || rc->rightmost);
    //     si->print();
    // 
    // }
    // 
    // free_ResultContig(rc);
    // free_ResultContig(crc);
    // free_CSList(root);
    // 
}
