#include <assert.h>
#include <errno.h>

#include "synmap.h"

Synmap * init_synmap(){
    Synmap* syn = (Synmap*)malloc(sizeof(Synmap));
    syn->genome = (Genome**)malloc(2 * sizeof(Genome*));
    return(syn);
}

void free_synmap(Synmap * synmap){
    if(synmap){
        free_genome(SG(synmap, 0));
        free_genome(SG(synmap, 1));
        free(synmap->genome);
        free(synmap);
    }
}

void print_synmap(Synmap * synmap){
    print_genome(SG(synmap, 0));
    print_genome(SG(synmap, 1));
}
