#include <assert.h>
#include <errno.h>

#include "global.h"
#include "synmap.h"
#include "genome.h"
#include "contig.h"
#include "block.h"

Synmap * init_synmap(){
    Synmap* syn = (Synmap*)malloc(sizeof(Synmap));
    syn->genome = (Genome**)malloc(2 * sizeof(Genome*));
    return(syn);
}

void free_synmap(Synmap * synmap){
    if(synmap){
        free_genome(synmap->genome[0]);
        free_genome(synmap->genome[1]);
        free(synmap->genome);
        free(synmap);
    }
}

void print_synmap(Synmap * synmap){
    print_genome(synmap->genome[0]);
    print_genome(synmap->genome[1]);
}

void check_args(int line_no, int nargs, int correct_nargs){
    if(nargs < correct_nargs){
        fprintf(stderr, "Either too few fields on input line %d,", line_no);
        fprintf(stderr, "or they are of the wrong type\n");
        exit(EXIT_FAILURE);
    }
    else if(nargs > correct_nargs){
        fprintf(stderr, "Too many fields on input line %d\n", line_no);
        exit(EXIT_FAILURE);
    }
}

Synmap * load_synmap(FILE * synfile){
    assert(synfile != NULL);

    Synmap * synmap = init_synmap();

    int line_no = 0;
    int unloaded_blocks = 0;
    int status;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char seqid[128];
    uint ncontigs, nblocks, id;

    // A dummy variable searched for at the end of each sscanf format sequence.
    // It will be matched only if there are too many arguments.
    char dummy;
    for(int i = 0; i < 2; i++){
        while ((read = getline(&line, &len, synfile)) != EOF) {
            line_no++;
            if(line[0] == '>'){
                status = sscanf(line, "> %s %u %c", seqid, &ncontigs, &dummy);
                check_args(line_no, status, 2);
                synmap->genome[i] = init_genome(seqid, ncontigs);
            }
            else if(line[0] == '@'){
                break;
            }
            else if(line[0] == '$') {
                status = sscanf(line, "$ %u %u %s %c\n", &id, &nblocks, seqid, &dummy);
                check_args(line_no, status, 3);
                synmap->genome[i]->contig[id] = init_contig(seqid, nblocks);
                unloaded_blocks += nblocks;
            }
            else {
                fprintf(stderr, "Incorrect file format, line %d\n", line_no);
                exit(EXIT_FAILURE);
            }
        }
    }

    line_no = 0;
    uint qcon_id, qblk_id, qstart, qstop;
    uint tcon_id, tblk_id, tstart, tstop;
    uint link_id;
    Block * new_block;
    while ((read = getline(&line, &len, synfile)) != EOF) {
        line_no++;
        if(line[0] != '$')
            continue;
        unloaded_blocks -= 2;
        status = sscanf(line, "$ %u %u %u %u %u %u %u %u %u %c\n",
                     &qcon_id, &qblk_id, &qstart, &qstop,
                     &tcon_id, &tblk_id, &tstart, &tstop,
                     &link_id, &dummy);
        check_args(line_no, status, 9);
        if(qstart > qstop || tstart > tstop){
            fprintf(stderr, "start must be less than stop on line %d\n", line_no);
            exit(EXIT_FAILURE);
        }
        // don't exceed the specified number of Contig in Genome
        if(qcon_id >= synmap->genome[0]->size || tcon_id >= synmap->genome[1]->size){
            fprintf(stderr, "too few contigs specified\n");
            exit(EXIT_FAILURE);
        }
        // don't exceed the specified size of the Contig block arrays
        if(qblk_id >= synmap->genome[0]->contig[qcon_id]->size ||
           tblk_id >= synmap->genome[1]->contig[tcon_id]->size){
            fprintf(stderr, "too few blocks specified\n");
            exit(EXIT_FAILURE);
        }

        new_block = init_block(qstart, qstop, tcon_id, tblk_id, link_id);
        synmap->genome[0]->contig[qcon_id]->block[qblk_id] = new_block;

        new_block = init_block(tstart, tstop, qcon_id, qblk_id, link_id);
        synmap->genome[1]->contig[tcon_id]->block[tblk_id] = new_block;
    }
    free(line);

    /* each block specified has been loaded */
    if(unloaded_blocks != 0){
        fprintf(stderr, "Unequal number of specified and actual blocks\n");
        exit(EXIT_FAILURE);
    }

    return(synmap);
}
