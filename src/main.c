#include <stdlib.h>

#include "global.h"
#include "ui.h"
#include "synmap.h"


int main(int argc, char * argv[]){

    Synmap * syn = NULL;

    // ------------------------------------------------------------------------
    // Prep input 
    // ------------------------------------------------------------------------
    Arguments args = parse_command(argc, argv);
    if(args.db_filename){
        char cmd[256];
        sprintf(cmd,
                "./util/prepare-data.sh -a %s -b %s -i %s -d %s",
                "susan", "philus", args.db_filename, "db");
        system(cmd);
    }
    if(args.synfile){
        syn = load_synmap(args.synfile);
    }

    // ------------------------------------------------------------------------
    // Do stuff 
    // ------------------------------------------------------------------------
    if(args.intfile && args.synfile){
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


    // ------------------------------------------------------------------------
    // Clean up
    // ------------------------------------------------------------------------
    free_synmap(syn);
    close_Arguments(args);

    return(EXIT_SUCCESS);
}
