#ifndef __INTERVAL_RESULT_H__
#define __INTERVAL_RESULT_H__

#include "global.h"

#include <vector>

// Forward declaration
template <class T> class IntervalTree; 

/** A container for results from searches for overlaps
 *
 * bool inbetween - TRUE if the query overlaps nothing
 */
template <class T>
class IntervalResult
{
public:
    std::vector<T*> iv;
    IntervalTree<T>* tree = nullptr;
    bool inbetween        = false;
    bool leftmost         = false;
    bool rightmost        = false;

    IntervalResult(){ };
    ~IntervalResult(){ }
};

#endif
