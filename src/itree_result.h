#ifndef __ITREE_RESULT_H__
#define __ITREE_RESULT_H__

#include <vector>

#include "interval.h"

/** A container for results from searches for overlaps
 *
 * bool inbetween - TRUE if the query overlaps nothing
 */
class IntervalResult
{
public:
    std::vector<Interval*> iv;
    bool inbetween;
    bool leftmost;
    bool rightmost;

    IntervalResult();

    void print();
};

#endif
