#ifndef __MANY_CONTIGUOUS_SETS_H__
#define __MANY_CONTIGUOUS_SETS_H__

#include "global.h"
#include "contiguous_set.h"
#include "interval_set.hpp"

/** A containter for ContiguousSets */
class ManyContiguousSets : public IntervalSet<ContiguousSet>
{
public:
    ManyContiguousSets();
    ~ManyContiguousSets();
    void link_contiguous_blocks(
        Block*  front,
        long    k,
        size_t& setid
    );

    /** Build a contiguous set from the homologous set
     *
     * @param first - The first block in the homologous set
     *
     */
    void add_from_homolog(ContiguousSet* first);
};

#endif
