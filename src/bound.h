#ifndef __BOUND_H__
#define __BOUND_H__

#include "interval.h"

class Bound : public Interval<Bound>
{
public:
    Bound();
    Bound(long start, long stop);
};

#endif
