#ifndef __RESULT_H__
#define __RESULT_H__

#include "global.h"

/**
 * Store results from synteny search
 *
 * This structure currently only contains the most basic return values.
 * Later I my expand this to cover whatever other metrics I gather.
 */
typedef struct {
    uint qseqid;
    uint qstart;
    uint qstop;
    uint tseqid;
    uint tstart;
    uint tstop;
    bool mismatch;
} Result;

#endif
