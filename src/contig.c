#include <string.h>

#include "global.h"
#include "contig.h"
#include "block.h"
#include "result.h"

#include "itree/itree.h"

/** Allocate memory for a contig and set each field.
 *
 * The itree field, which can hold an interval tree for (log(n)+m) overlap
 * searching, is initialized to NULL. Building this tree is expensive (n
 * log(n)) so will not be done unless needed.
 *
 * Memory is allocated for pointers to *size* blocks, but they are not
 * initialized.
 *
 * @param name name of this Contig (e.g. "Chr1")
 * @param size number of Block structs this Contig will hold
 *
 * @return pointer to a new Contig
 *
 * */
Contig * init_contig(char * name, size_t size){
    Contig* con = (Contig*)malloc(sizeof(Contig));
    con->name = strdup(name);
    con->size = size;
    con->itree = NULL;
    con->block = (Block**)malloc(size * sizeof(Block*));
    return con;
}

/** Recursively free all memory.
 *
 * This functions calls free_block on each Block in its block field.
 *
 * If the Contig has an IntervalTree defined, it will free it with free_interval_tree.
 *
 * @param contig pointer to a contig, may be NULL
 * */
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

/** Recursively print contig. */
void print_contig(Contig * contig){
    printf("%lu\t%s\n", contig->size, contig->name);
    for(int i = 0; i < contig->size; i++){
        print_block(contig->block[i]);
    }
}

uint anchor(uint x, Contig * contig){
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

Contig * get_overlapping(uint a, uint b, Contig * con){
    // STUB
    return con;
}

Contig * get_flanks(uint a, uint b, Contig * con, uint nup, uint ndown){
    // STUB
    return con;
}

uint count_overlaps(uint a, uint b, Contig * con){
    // STUB
    return a;
}

Result map(uint a, uint b, Contig * con){
    // STUB
    Result mr;
    return mr;
}
