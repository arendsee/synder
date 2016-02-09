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

// local function
IA * ia_from_blocks(Contig * con){
    IA * ia = init_set_ia(con->size);
    for(int i = 0; i < con->size; i++){
        ia->v[i].start = con->block[i]->start;
        ia->v[i].stop = con->block[i]->stop;
        ia->v[i].link = (void*)con->block[i];
    }
    return ia;
}

Contig * get_overlapping(Contig * con, uint a, uint b){
    if(!con->itree)
        con->itree = build_tree(ia_from_blocks(con));
    Interval * inv = init_interval(a, b); 
    IA * ia = get_interval_overlaps(inv, con->itree);
    Contig * newcon = init_contig(con->name, ia->size);
    for(int i = 0; i < ia->size; i++){
        newcon->block[i] = (Block*)(ia->v[i].link);
    }
    free_ia(ia);
    free(inv);
    return newcon;
}

Contig * get_flanks(Contig * con, uint a, uint b, uint nup, uint ndown){
    // STUB
    return con;
}

uint count_overlaps(Contig * con, uint a, uint b){
    if(!con->itree)
        con->itree = build_tree(ia_from_blocks(con));
    Interval * inv = init_interval(a, b); 
    uint count = count_interval_overlaps(inv, con->itree);
    free(inv);
    return count;
}

Result map(Contig * con, uint a, uint b){
    // STUB
    Result mr;
    return mr;
}
