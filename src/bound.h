#ifndef __BOUND_H__
#define __BOUND_H__

#include "interval.h"

class Bound : Interval<Bound>
{
    Bound(long start, long stop);
};

#endif
