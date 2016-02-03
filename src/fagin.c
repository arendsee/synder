#include "ui.h"
#include "syndb.h"
#include "ftypes.h"

int main(int argc, char * argv[]){

    // Prep input
    Arguments args = parse_command(argc, argv);
    Synmap * syn = load_synmap(args.synfile);

    // Do stuff
    print_synmap(syn);

    // Cleanup
    free_synmap(syn);
    close_Arguments(args);

    return(EXIT_SUCCESS);
}
