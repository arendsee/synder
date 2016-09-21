#include "many_search_intervals.h"

ManySearchIntervals::~ManySearchIntervals{
    for(auto &i : inv){
        delete i;
    }
}

void ManySearchIntervals::find_search_intervals(const Feature& t_feat){
    // Does starget strand == '-'?
    bool inverted;

    std::set<ContiguousSet*> csets;

    std::unique_ptr< IntervalResult<Block> > rc = cset->get_region(t_feat, true);

    // get list of highest and lowest members of each contiguous set
    for (auto &q : rc->iv) {
        csets.insert(q->cset);
    }

    // TODO what am I doing here?
    // Merge all this crap into the SearchInterval class
    std::unique_ptr< IntervalResult<ContiguousSet> > crc = cset->get_region(t_feat, false);
    if(! (crc->inbetween || crc->leftmost || crc->rightmost) ) {
        for (auto &q : crc->iv) {
            csets.insert(q);
        }
    }

    // Iterate through each contiguous set, for each find the search interval
    // For each contiguous set, there is exactly one search interval, or into a new SearchIntervalSet class
    for(auto &cset : csets) {
        bool inbetween = rc->inbetween || rc->leftmost || rc->rightmost
        inv.push_back( new SearchInterval(cset->ends, &t_feat, inbetween) );
    }

    for(auto &i : inv){
        i.print();
    }
}
