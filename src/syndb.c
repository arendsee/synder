#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "syndb.h"
#include "itree/itree.h"

// Local initiator funcions
Block * init_block(uint, uint, uint, uint, uint);
Contig * init_contig(char *, size_t);
Genome * init_genome(char *, size_t);
Synmap * init_synmap();

// Local memory freeing functions
void free_block(Block *);
void free_contig(Contig *);
void free_genome(Genome *);


/** Recursively free all memory allocated to the synteny map.
 * 
 * Calls free_genome on both its Genome children.
 *
 * @param pointer to Synmap struct
 *
 * */
void free_synmap(Synmap * synmap){
    free_genome(synmap->genome[0]);
    free_genome(synmap->genome[1]);
    free(synmap->genome);
    free(synmap);
}

/** Recursively print a synteny map. */
void print_synmap(Synmap * synmap){
    print_genome(synmap->genome[0]);
    print_genome(synmap->genome[1]);
}

/** Recursively print a genome. */
void print_genome(Genome * genome){
    printf(">\t%s\t%lu\n", genome->name, genome->size);
    for(int i = 0; i < genome->size; i++){
        print_contig(genome->contig[i]);
    }
}

/** Recursively print contig. */
void print_contig(Contig * contig){
    printf("%lu\t%s\n", contig->size, contig->name);
    for(int i = 0; i < contig->size; i++){
        print_block(contig->block[i]);
    }
}

/** Print all fields in this block (TAB-delimited). */
void print_block(Block * block){
    printf("%u\t%u\t%u\t%u\t%u\n", 
        block->start,
        block->stop,
        block->oseqid,
        block->oblkid,
        block->linkid);
}

/** Build synteny tree from specially formatted file.
 *
 * @warning This function is VERY picky about input. It expects input to be
 * formatted exactly as util/prepare-data.sh produces. You must not feed this
 * function raw synteny files. I currently have no input checks.
 *
 * @param synfile specially formatted synteny file
 *
 * @return pointer to a complete Synmap object
 */

Synmap * load_synmap(FILE * synfile){
    assert(synfile != NULL);

    Synmap * synmap = init_synmap();

    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char seqid[128];
    uint ncontigs, nblocks, id;
    for(int i = 0; i < 2; i++){
        while ((read = getline(&line, &len, synfile)) != EOF) {
            if(line[0] == '>'){
                sscanf(line, "> %s %u", seqid, &ncontigs);
                synmap->genome[i] = init_genome(seqid, ncontigs);
            }
            else if(line[0] == '@'){
                break;
            } else {
                sscanf(line, "%u %u %s", &id, &nblocks, seqid);
                synmap->genome[i]->contig[id] = init_contig(seqid, nblocks);
            }
        }
    }
    free(line);

    uint qcon_id, qblk_id, qstart, qstop;
    uint tcon_id, tblk_id, tstart, tstop;
    uint link_id;
    Block * new_block;
    while ((fscanf(synfile, "%u %u %u %u %u %u %u %u %u",
                   &qcon_id, &qblk_id, &qstart, &qstop,
                   &tcon_id, &tblk_id, &tstart, &tstop,
                   &link_id
            ) != EOF)) {
        new_block = init_block(qstart, qstop, tcon_id, tblk_id, link_id);
        synmap->genome[0]->contig[qcon_id]->block[qblk_id] = new_block;

        new_block = init_block(tstart, tstop, qcon_id, qblk_id, link_id);
        synmap->genome[1]->contig[tcon_id]->block[tblk_id] = new_block;
    }

    return(synmap);
}



/******************************************************************************
 * Local functions
 *****************************************************************************/

/** Allocate memory for a block and set each field.
 *
 * @param start  query start position
 * @param stop   query stop position
 * @param oseqid index of target Contig
 * @param oblkid index of matching Block on target Contig
 * @param linkid index of this Block pair in metadata array
 *
 * @return pointer to new Block
 *
 * */
Block * init_block(uint start, uint stop,
                   uint oseqid, uint oblkid,
                   uint linkid){
    Block * block = (Block*)malloc(sizeof(Block));
    block->start  = start;
    block->stop   = stop;
    block->oseqid = oseqid;
    block->oblkid = oblkid;
    block->linkid = linkid;
    return(block);
}
/** Free memory allocated to this block.
 *
 * @param block pointer to a Block, may be NULL
 * */
void free_block(Block * block){
    if(block){
        free(block);
    }
}


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


/** Allocate memory for Genome *name* of size *size*.
 *
 * @param name genome name (e.g. "Arabidopsis_thaliana")
 * @param size number of child Contig structs (e.g. chromosomes or scaffolds) 
 *
 * @return pointer to new Genome struct 
 * */
Genome * init_genome(char * name, size_t size){
    Genome * gen = (Genome*)malloc(sizeof(Genome));
    gen->name = strdup(name);
    gen->size = size;
    gen->contig = (Contig**)malloc(size * sizeof(Contig*));
    return(gen);
}
/**
 * Recursively tree all memory
 *
 * For each Contig in the contig field, calls free_contig.
 *
 * @param pointer to a Genome struct
 */
void free_genome(Genome * genome){
    if(genome){
        for(int i = 0; i < genome->size; i++){
            free_contig(genome->contig[i]);
        }
        free(genome->contig);
        free(genome->name);
        free(genome);
    }
}


/**
 * Allocate memory for a new Synmap.
 *
 * This struct holds the pair of genomes that wil be compared.
 *
 * @return pointer to the new Synmap
 */
Synmap * init_synmap(){
    Synmap* syn = (Synmap*)malloc(sizeof(Synmap));
    syn->genome = (Genome**)malloc(2 * sizeof(Genome*));
    return(syn);
}
