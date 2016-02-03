#include "ftypes.h"
#include "ui.h"
#include "syndb.h"
#include "interval.h"

int main(int argc, char * argv[]){

    // Prep input
    Arguments args = parse_command(argc, argv);
    Synmap * syn = load_synmap(args.synfile);

    // Do stuff
    // print_synmap(syn);
    uint x = 15000;
    uint anc = anchor(x, syn->genome[0]->contig[0]);
    printf("%u\n", anc);

    // Cleanup
    free_synmap(syn);
    close_Arguments(args);

    return(EXIT_SUCCESS);
}
