#include <string.h>

#include "global.h"
#include "contig.h"
#include "block.h"
#include "result.h"

#include "itree/itree.h"

Contig * init_contig(char * name, size_t size){
    Contig* con = (Contig*)malloc(sizeof(Contig));
    con->name = strdup(name);
    con->size = size;
    con->itree = NULL;
    con->block = (Block**)malloc(size * sizeof(Block*));
    return con;
}

void free_contig(Contig * contig){
    if(contig){
        if(contig->itree)
            free_interval_tree(contig->itree);
        for(int i = 0; i < contig->size; i++){
            free_block(contig->block[i]);
        }
        free(contig->block);
        free(contig->name);
        free(contig);
    }
}

void print_contig(Contig * contig){
    printf("%lu\t%s\n", contig->size, contig->name);
    for(int i = 0; i < contig->size; i++){
        print_block(contig->block[i]);
    }
}

uint anchor(Contig * contig, uint x){
    Block ** blks = contig->block;
    uint N = contig->size;
    uint lo = 0;
    uint hi = N - 1;
    uint i = hi / 2;
    while(true){
        if(x >= blks[i]->start){
            if(i == (N-1) || x < blks[i+1]->start)
                return i;
            lo = i;
            i = (i + hi) / 2 + 1;
        }
        else{
            if(i == 0 || x > blks[i-1]->start)
                return i == 0 ? i : i-1;
            hi = i;
            i = lo + (i - lo) / 2;
        }
    }
    return i;
}

Contig * get_overlapping(Contig * con, uint a, uint b){
    // STUB
    return con;
}

Contig * get_flanks(Contig * con, uint a, uint b, uint nup, uint ndown){
    // STUB
    return con;
}

uint count_overlaps(Contig * con, uint a, uint b){
    // STUB
    return a;
}

Result map(Contig * con, uint a, uint b){
    // STUB
    Result mr;
    return mr;
}
