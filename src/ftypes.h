#ifndef __FTYPES_H__
#define __FTYPES_H__

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define uint unsigned int

typedef struct {
    uint start;
    uint stop;
    uint oseqid;
    uint oblkid;
    uint linkid;
} Block;

typedef struct {
    char * name;
    size_t size;
    Block ** block;
} Contig;

typedef struct {
    char * name;
    size_t size;
    Contig ** contig;
} Genome;

typedef struct {
    Genome ** genome;
} Synmap;

#endif
