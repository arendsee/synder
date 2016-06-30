#ifndef __SYNMAP_H__
#define __SYNMAP_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#define TOLC(s)	for(char *p =(s);*p;++p) *p=*p>0x40&&*p<0x5b?*p|0x60:*p
#define LEVMIN(x,y,z) ((x)<(y)?((x)<(z)?(z):(z)):((y)<(z)?(y):(z)))

typedef struct  DictNode{
    char * word;
    uint32_t arr_loc;
	bool original;
    struct DictNode * next;
} DictNode;

int32_t lev_dist(char *s1, char *s2);

void convert_seqname(FILE* synfile, FILE* intfile);

#endif
