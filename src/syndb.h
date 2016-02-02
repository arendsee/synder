#ifndef __SYNDB_H__
#define __SYNDB_H__

typedef struct {
    unsigned int start;
    unsigned int stop;
    unsigned int oseqid;
    unsigned int oblkid;
    unsigned int linkid;
} interval;

typedef struct {
    size_t size;
    block * blocks;
} contig

typedef struct {
    size_t size;
    contig * contigs;
} contig_set 

typedef struct {
    char * name_a;
    size_t size_a;
    contig * a;
    char * name_b;
    size_t size_b;
    contig * b;
} synteny_pair;

#endif
