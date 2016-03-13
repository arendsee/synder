#ifndef __CONTIG_H__
#define __CONTIG_H__

struct {
    size_t size;
    struct S_Interval ** intervals;
} S_Contig

struct S_Contig init_s_contig();

void free_s_contig(struct S_Contig * contig);

void print_s_contig(struct S_Contig * contig);

#endif
