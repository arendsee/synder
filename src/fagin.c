#include <stdio.h>
#include <stdlib.h>

#include "ui.h"
#include "syndb.h"

int main(int argc, char * argv[]){

    Arguments args = parse_command(argc, argv);

    Synmap * syn = load_synmap(args.synfile);

    // Closes synteny map file if open
    close_Arguments(args);

    free_synmap(syn);

    return(EXIT_SUCCESS);
}
