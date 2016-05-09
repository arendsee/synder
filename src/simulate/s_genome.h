#ifndef __CONTIG_H__
#define __CONTIG_H__

struct {
    size_t size;
    struct S_Contig ** contigs;
} S_Genome

struct S_Genome init_s_genome();

void print_s_genome(struct S_Genome * genome);

void free_s_genome(struct S_Genome * genome);

#endif
