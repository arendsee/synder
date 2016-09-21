#ifndef __SEARCH_INTERVAL_H__
#define __SEARCH_INTERVAL_H__

#include "interval.h"
#include "feature.h"

enum Flag {
    ANCHORED = 0, // bound in inside a syntenic interval
    BOUND    = 1, // bound is between members of a contiguous set
    UNBOUND  = 2, // bound is between contiguous sets
    BEYOND   = 3, // bound is further out than any syntenic interval
    EXTREME  = 4  // bound is BEYOND and is tangent to the contig bound
};

class SearchInterval : Interval<SearchInterval>
{
private:

    Feature*             m_feat      = nullptr;
    bool                 m_inbetween = false;
    std::array<Block*,2> m_bnds      = {     nullptr, nullptr };
    double               m_score     = 0;
    std::array<int,2>    m_flag      = {404, 404};
    bool                 m_inverted  = false;

    void set_bound(Direction d); 
    double flank_area(long near, long far, double k)
    double calculate_score(Bound& bound, Block* blk)

public:
    SearchInterval(const Feature& feat);
    ~SearchInterval();

    void build_search_interval(
        ContiguousSet* t_cset,
        const Feature& t_feat,
        bool t_inbetween
    );

    void print();

}

#endif
