#ifndef __INTERVAL_H__
#define __INTERVAL_H__

struct {
    uint start;
    uint stop;
    bool strand;
    size_t nlinks;
    struct S_Interval ** links;
} S_Interval;


struct S_Interval init_s_interval();

void free_s_interval(struct S_Interval * interval);

void print_s_interval(struct S_Interval * interval);

#endif
