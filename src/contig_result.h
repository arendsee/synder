#ifndef __CONTIG_RESULT_H__
#define __CONTIG_RESULT_H__

#include "global.h"
#include "itree.h"


/** A wrapper for Contig that includes IntervalResult flags
 *
 * Fields
 *  - contig - a pointer to the original Contig
 *  - size   - number of elements in block
 *  - block  - pointer to selection of Block structs in Contig
 *  - flags
 *    1. inbetween - query is between (not overlapping) two search intervals
 *    2. leftmost  - query is further left than any interval
 *    3. rightmost - query is further right than any interval
 */
typedef struct ResultContig
{
    Contig* contig;
    size_t size;
    Block** block;
    ContiguousSet** cset;
    bool inbetween;
    bool leftmost;
    bool rightmost;
} ResultContig;

ResultContig* init_ResultContig(Contig*, IntervalResult*, bool is_cset);

void free_ResultContig(ResultContig*);

#endif
