#include "interval.h"

Interval::Interval(
    long new_start,
    long new_stop
)
    : start(new_start), stop(new_stop)
{
    link = NULL;
}

void Interval::print()
{
    printf("%zu %zu\n", start, stop);
}

bool Interval::cmp_start(Interval * a, Interval * b)
{
    return ( a->start < b->start );
}

bool Interval::cmp_stop(Interval * a, Interval * b)
{
    return ( a->stop < b->stop );
}

Pos Interval::overlap(long a)
{
    if (a < start)
    {
        return lo;
    }
    else if (a > stop)
    {
        return hi;
    }
    else
    {
        return in;
    }
}

Pos Interval::overlap(Interval * b)
{
    if (stop < b->start)
    {
        return lo;
    }
    else if (start > b->stop)
    {
        return hi;
    }
    else
    {
        return in;
    }
}
