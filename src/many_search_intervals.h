#ifndef __MANY_SEARCH_INTERVALS_H__
#define __MANY_SEARCH_INTERVALS_H__

#include <map>

class ManySearchIntervals : public IntervalSet<SearchInterval>
{
public:
    ManySearchIntervals();
    ~ManySearchIntervals();
};

#endif
