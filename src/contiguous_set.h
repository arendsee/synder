#ifndef __CONTIGUOUS_SET_H__
#define __CONTIGUOUS_SET_H__

#include "global.h"
#include "block.h"


ContiguousSet * init_ContiguousSet(Block * blk);

bool add_block_to_ContiguousSet(ContiguousSet * cset, Block * blk, long k);

/** Determine whether two blocks are contiguous
 * 
 * Each block consists of query and target intervals
 * 
 * Let these 4 resulting intervals be aq, at, bq and bt
 * 
 * Each interval is defined by 3 variables:
 *   1. s - target strand
 *   2. c - target chromosome/scaffold name
 *   3. g - non-overlapping group id
 * 
 * I will identify each of these variables by appending [scg] to the
 * interval name, e.g. ats or bqg.
 * 
 * The blocks are contiguous if and only if all of the following are true
 *   1. aqs == bqs
 *   2. ats == bts
 *   3. aqc == bqc
 *   4. atc == btc
 * 
 *   #1 will always be true, since strand is relative to query.
 */
bool are_contiguous(
    Block * blk_a,
    Block * blk_b,
    long k
);

bool strictly_forbidden(
    Block * a,
    Block * b,
    long k
);

#endif
