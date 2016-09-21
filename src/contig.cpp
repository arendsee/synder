#include "contig.h"

Contig::Contig(const char* t_genome_name, const char* t_contig_name, long t_length)
:
    feat(t_genome_name, 0, t_length, t_contig_name, t_length)
{ }

void Contig::set_length(long t_length)
{
    feat.parent_length = t_length;
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
    // for (; c != nullptr; c = c->next) {
    //     fprintf(stderr, "  -- ");
    //     print_ContiguousSet(c);
    // }
    // 
    // if (print_blocks) {
    //     int d = forward ? 1 : 3; // next by start or next by stop
    //     Block* blk = cor[d - 1]; // prev by start or prev by stop
    //     for (; blk != nullptr; blk = blk->cor[d]) {
    //         print_Block(blk);
    //     }
    // }
}

void Contig::count(Feature& t_feat)
{
    long count = block.count_overlaps(&t_feat);
    printf("%s\t%zu\n", t_feat.name.c_str(), count);
}

void Contig::map(Feature& t_feat)
{
    auto rc = block.get_region(t_feat, true);
    bool missing = rc->inbetween || rc->leftmost || rc->rightmost;

    for (auto &qblk : rc->iv) {
        if (qblk != nullptr) {
            printf("%s %s %zu %zu %d\n",
                   t_feat.name.c_str(),
                   t_feat.parent_name.c_str(),
                   qblk->over->pos[0] + Offsets::out_start,
                   qblk->over->pos[1] + Offsets::out_stop,
                   missing
                  );
        }
    }
    delete rc;
}

void Contig::find_search_intervals(Feature& t_feat){
    
}
