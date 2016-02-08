#include <assert.h>

#include "global.h"
#include "synmap.h"
#include "genome.h"
#include "contig.h"
#include "block.h"

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
