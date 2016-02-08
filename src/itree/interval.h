#ifndef __INTERVAL_H__
#define __INTERVAL_H__

/* the eponymous structure */
typedef struct {
    unsigned int start;
    unsigned int stop;
} Interval;

/* position of Interval or point A relative to B */
typedef enum {lo, in, hi} Pos;

/* compare intervals by stop */
int cmp_stop(const void *, const void *);

/* compare intervals by start */
int cmp_start(const void *, const void *);

/* find position of point A relative to interval B (see Pos) */
Pos point_overlap(unsigned int A, Interval B);

/* find position of interval A relative to interval B (see Pos) */
Pos overlap(Interval, Interval);

#endif
