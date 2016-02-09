#include "global.h"
#include "ui.h"
#include "synmap.h"

int main(int argc, char * argv[]){

    // ------------------------------------------------------------------------
    // Prep input 
    // ------------------------------------------------------------------------
    Arguments args = parse_command(argc, argv);
    Synmap * syn = load_synmap(args.synfile);


    // ------------------------------------------------------------------------
    // Do stuff 
    // ------------------------------------------------------------------------

    //print_synmap(syn);
    print_args(args);

    uint count = count_overlaps(syn->genome[0]->contig[args.chr], args.start, args.stop);
    printf("%u\n", count);


    // ------------------------------------------------------------------------
    // Clean up
    // ------------------------------------------------------------------------
    free_synmap(syn);
    close_Arguments(args);

    return(EXIT_SUCCESS);
}
