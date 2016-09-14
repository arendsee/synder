#include "itree_result.h"

IntervalResult::IntervalResult()
{
    inbetween = false;
    leftmost  = false;
    rightmost = false;
}

void IntervalResult::print()
{
    printf("inbetween=%i leftmost=%i rightmost=%i\n", inbetween, leftmost, rightmost);
}
