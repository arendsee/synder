#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "syndb.h"

Block * init_block(uint start, uint stop,
                   uint oseqid, uint oblkid,
                   uint linkid){
    Block * block = (Block*)malloc(sizeof(Block));
    block->start = start;
    block->stop = stop;
    block->oseqid = oseqid;
    block->oblkid = oblkid;
    block->linkid = linkid;
    return(block);
}

void free_block(Block * block){
    if(block){
        free(block);
    }
}

Contig * init_contig(char * name, size_t size){
    Contig* con = (Contig*)malloc(sizeof(Contig));
    con->name = strdup(name);
    con->size = size;
    con->block = (Block**)malloc(size * sizeof(Block*));
    return con;
}

void free_contig(Contig * contig){
    if(contig){
        for(int i = 0; i < contig->size; i++){
            free_block(contig->block[i]);
        }
        free(contig->block);
        free(contig->name);
        free(contig);
    }
}

Genome * init_genome(char * name, size_t size){
    Genome * gen = (Genome*)malloc(sizeof(Genome));
    gen->name = strdup(name);
    gen->size = size;
    gen->contig = (Contig**)malloc(size * sizeof(Contig*));
    return(gen);
}

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

Synmap * init_synmap(){
    Synmap* syn = (Synmap*)malloc(sizeof(Synmap));
    syn->genome = (Genome**)malloc(2 * sizeof(Genome*));
    return(syn);
}

void free_synmap(Synmap * synmap){
    free_genome(synmap->genome[0]);
    free_genome(synmap->genome[1]);
    free(synmap->genome);
    free(synmap);
}

Synmap * load_synmap(FILE * synfile){
    Synmap * synmap = init_synmap();

    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char seqid[128];
    uint ncontigs;
    uint nblocks;
    uint id;
    assert(synfile != NULL);
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

    uint qcon_id;
    uint qblk_id;
    uint qstart;
    uint qstop;
    uint tcon_id;
    uint tblk_id;
    uint tstart;
    uint tstop;
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
