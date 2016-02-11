#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "ui.h"
#include "synmap.h"
#include "analysis.h"

int main(int argc, char * argv[]){

    Synmap * syn = NULL;

    // ------------------------------------------------------------------------
    // Prep input 
    // ------------------------------------------------------------------------

    Arguments args = parse_command(argc, argv);

    // ------------------------------------------------------------------------
    // Do stuff 
    // ------------------------------------------------------------------------

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

    if(args.synfile){
        syn = load_synmap(args.synfile);
    }

    if(!syn){
        printf("Nothing to do ...\n");
        print_help();
    }

    if(args.intfile){
        if(strcmp(args.cmd, "count") == 0){
            analysis_count(syn, args.intfile);
        }
        else if(strcmp(args.cmd, "map") == 0){
            analysis_map(syn, args.intfile);
        }
        else{
            printf("Command '%s' not recognized\n", args.cmd);
            print_help();
        }
    }


    // ------------------------------------------------------------------------
    // Clean up
    // ------------------------------------------------------------------------
    free_synmap(syn);
    close_Arguments(args);

    return(EXIT_SUCCESS);
}
