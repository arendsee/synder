#include <stdio.h>
#include <stdlib.h>

#include "ui.h"

int main(int argc, char * argv[]){

    Arguments args = parse_command(argc, argv);

    print_args(args);

    close_Arguments(args);

    return(EXIT_SUCCESS);
}
