#include "contig.h"

Contig::Contig() { }

Contig::Contig(const char* t_genome_name, const char* t_contig_name, long t_length)
    :
    feat(t_genome_name, 0, t_length, t_contig_name, t_length)
{ }

Contig::~Contig() { }

void Contig::set_length(long t_length)
{
    feat.parent_length = t_length;
}

void Contig::count(Feature& t_feat, CountType& df)
{
    long count = block.count_overlaps(&t_feat);

    df.add_row(t_feat.name, count);
}

void Contig::map(Feature& t_feat, MapType& df)
{
    auto rc = block.get_region(t_feat, true);
    bool missing = rc->inbetween || rc->leftmost || rc->rightmost;

    for (auto &qblk : rc->iv) {
        if (qblk != nullptr) {
            df.add_row(
                t_feat.name,
                qblk->parent->parent_name,
                qblk->pos[0] + Offsets::out_start,
                qblk->pos[1] + Offsets::out_stop,
                t_feat.parent_name,
                qblk->over->pos[0] + Offsets::out_start,
                qblk->over->pos[1] + Offsets::out_stop,
                qblk->over->strand,
                missing
            );
        }
    }
    delete rc;
}

std::vector<SearchInterval> Contig::list_search_intervals(Feature& t_feat, double r)
{
    // TODO -- need to move this back up to Contig

    std::set<ContiguousSet*> csets;

    auto rc = block.get_region(t_feat, true);

    // get list of highest and lowest members of each contiguous set
    for (auto &q : rc->iv) {
        csets.insert(q->cset);
    }

    // TODO what am I doing here?
    // Merge all this crap into the SearchInterval class
    auto crc = cset.get_region(t_feat, false);
    if(! (crc->inbetween || crc->leftmost || crc->rightmost) ) {
        for (auto &q : crc->iv) {
            csets.insert(q);
        }
    }

    // Iterate through each contiguous set, for each find the search interval
    // For each contiguous set, there is exactly one search interval, or into a new SearchIntervalSet class
    bool inbetween = rc->inbetween || rc->leftmost || rc->rightmost;
    std::vector<SearchInterval> si;
    for(auto &c : csets) {
        si.push_back( SearchInterval(c->ends, &t_feat, inbetween, r) );
    }

    delete rc;
    delete crc;

    return si;
}

void Contig::find_search_intervals(Feature& t_feat, double r, SIType& stype)
{
    std::vector<SearchInterval> si = list_search_intervals(t_feat, r);
    for(auto &s : si) {
        s.add_row(stype);
    }
}
