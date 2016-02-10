#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "ui.h"
#include "synmap.h"


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
    
    if(!syn){
        printf("Nothing to do ...\n");
        print_help();
    }
    else if(strcmp(args.cmd, "count") == 0 && args.intfile){
        uint count = 0;
        char seqname[128];
        int chrid, start, stop;
        while ((fscanf(args.intfile,
                       "%d %*s %*s %d %d %*c %*c %*s %s\n",
                       &chrid, &start, &stop, seqname)) != EOF) {
            count = count_overlaps(syn->genome[0]->contig[chrid], start, stop);
            printf("%s\t%u\n", seqname, count);
        }
    }
    else{
        printf("No command found\n");
        print_help();
    }


    // ------------------------------------------------------------------------
    // Clean up
    // ------------------------------------------------------------------------
    free_synmap(syn);
    close_Arguments(args);

    return(EXIT_SUCCESS);
}
