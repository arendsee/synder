#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "ui.h"
#include "synmap.h"

/** \todo Move all the analysis crap from main.c to dedicated file */

void analysis_count(Synmap * syn, FILE * intfile){
    char seqname[128];
    uint count;
    int chrid, start, stop;
    while ((fscanf(intfile,
                   "%d %*s %*s %d %d %*c %*c %*s %s\n",
                   &chrid, &start, &stop, seqname)) != EOF)
    {
        count = count_overlaps(syn->genome[0]->contig[chrid], start, stop);
        printf("%s\t%u\n", seqname, count);
    }
}

void analysis_map(Synmap * syn, FILE * intfile){
    char seqname[128];
    int chrid, start, stop;
    Contig * contigs;
    Contig * tcon;
    Block * qblk;
    Block * tblk;
    while ((fscanf(intfile,
                   "%d %*s %*s %d %d %*c %*c %*s %s\n",
                   &chrid, &start, &stop, seqname)) != EOF)
    {
        contigs = get_overlapping(syn->genome[0]->contig[chrid], start, stop);
        for(int i = 0; i < contigs->size; i++){
            qblk = contigs->block[i];
            tcon = syn->genome[1]->contig[qblk->oseqid];
            tblk = tcon->block[qblk->oblkid];
            printf("%s %s %u %u\n", seqname, tcon->name, tblk->start, tblk->stop);
        }
    }
}

int main(int argc, char * argv[]){

    Synmap * syn = NULL;

    // ------------------------------------------------------------------------
    // Prep input 
    // ------------------------------------------------------------------------
    Arguments args = parse_command(argc, argv);

    // Build database and exit
    if(args.db_filename){
        if(!(args.pos[0] && args.pos[1] && args.pos[2]))
            print_help();
        char cmd[256];
        sprintf(cmd,
                "prepare-data.sh -a %s -b %s -i %s -d %s",
                args.pos[0], args.pos[1], args.db_filename, args.pos[2]);
        system(cmd);
        exit(EXIT_SUCCESS);
    }

    // Load synteny map from database
    if(args.synfile){
        syn = load_synmap(args.synfile);
    }

    // ------------------------------------------------------------------------
    // Do stuff 
    // ------------------------------------------------------------------------
    
    typedef enum {COUNT, MAP} Command;
    Command cmd; 
    if(!args.cmd){
        printf("No command found\n");
        print_help();
    }
    else if(strcmp(args.cmd, "count") == 0){
        cmd = COUNT; 
    }
    else if(strcmp(args.cmd, "map") == 0){
        cmd = MAP;
    }
    else{
        printf("Command '%s' not recognized\n", args.cmd);
        print_help();
    }
    
    if(!syn){
        printf("Nothing to do ...\n");
        print_help();
    }

    if(args.intfile){
        switch(cmd){
            case COUNT:
                analysis_count(syn, args.intfile);
                break;
            case MAP:
                analysis_map(syn, args.intfile);
                break;
            default:
                printf("Invalid command bypassed security\n");
                exit(EXIT_FAILURE);
        }

    }


    // ------------------------------------------------------------------------
    // Clean up
    // ------------------------------------------------------------------------
    free_synmap(syn);
    close_Arguments(args);

    return(EXIT_SUCCESS);
}
