#include "contig.h"

Contig::Contig(const char* t_genome_name, const char* t_contig_name, long t_length);
:
    block(), cset(), feat(t_genome_name, 0, t_length, t_contig_name, t_length)
{ }

void Contig::set_length(long t_length)
{
    length = feat.length = t_length;
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

void Contig::count(const Feature& feat)
{
    build_block_itree();
    long count = itree->count_overlaps(feat);
    printf("%s\t%zu\n", seqname, count);
}

void Contig::map(const Feature& feat)
{
    Block *qblk, *tblk;
    ContigResult* rc = get_region(feat, true);
    bool missing = rc->inbetween || rc->leftmost || rc->rightmost;

    for (size_t i = 0; i < rc->size; i++) {
        qblk = rc->block[i];
        if (qblk != NULL) {
            tblk = qblk->over;
            printf("%s %s %zu %zu %d\n",
                   feat.feature_name,
                   feat.contig_name,
                   tblk->pos[0] + Offsets::out_start,
                   tblk->pos[1] + Offsets::out_stop,
                   missing
                  );
        }
    }
    free(rc);
}

void Contig::find_search_intervals(const Feature& feat){
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
    // Output search interval
    SearchInterval si(this);

    // Get blocks overlapping the query
    rc = block->get_region(feat, true);

    // get list of highest and lowest members of each contiguous set
    cslist = init_empty_CSList();
    root = cslist;
    for(size_t i = 0; i < rc->size; i++) {
        add_blk_CSList(cslist, rc->block[i]);
    }

    // TODO what am I doing here?
    // Merge all this crap into the SearchInterval class
    crc = cset->get_region(feat, false);
    if(! (crc->inbetween || crc->leftmost || crc->rightmost) ) {
        for(size_t i = 0; i < crc->size; i++) {
            add_cset_CSList(cslist, crc->cset[i], feat);
        }
    }

    // Iterate through each contiguous set, for each find the search interval
    // For each contiguous set, there is exactly one search interval, or into a new SearchIntervalSet class
    for(; cslist != NULL; cslist = cslist->next) {

        blk_bounds[LO] = cslist->bound[LO];
        blk_bounds[HI] = cslist->bound[HI];

        // TODO fix the smell
        si = build_search_interval(feat, blk_bounds, seqname, rc->inbetween || rc->leftmost || rc->rightmost);
        si->print();

    }

    free_ResultContig(rc);
    free_ResultContig(crc);
    free_CSList(root);

}
