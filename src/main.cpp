#include "global.h"
#include "arguments.h"
#include "synmap.h"

// Global variables: 0/1 bases for input and output
int Offsets::in_start;
int Offsets::in_stop;
int Offsets::out_start;
int Offsets::out_stop;

int main(int argc, char *argv[])
{

    int exit_status = 0;
    Synmap* syn = NULL;

    // ------------------------------------------------------------------------
    // Prep input
    // ------------------------------------------------------------------------

    Arguments args = Arguments(argc, argv);

    Offsets::in_start  = args.offsets[0];
    Offsets::in_stop   = args.offsets[1];
    Offsets::out_start = args.offsets[2];
    Offsets::out_stop  = args.offsets[3];

    // ------------------------------------------------------------------------
    // Do stuff
    // ------------------------------------------------------------------------

    // If no file given by -i, use STDIN
    if (args.intfile == NULL)
        args.intfile = stdin;

    if (args.synfile) {
        syn = new Synmap(args);

        bool command_success = true;
        switch(args.cmd){
            case C_DEBUG:
                args.print();
                syn->print(true);
                break;
            case C_DUMP:
                syn->dump_blocks();
                break;
            default:
                command_success = syn->process_gff(args.intfile, args.cmd);
        }

        if(!command_success){
            fprintf(stderr, "Command failed\n");
            exit_status = 1;
        }

        delete syn;
    } else {
        printf("Nothing to do ...\n");
        args.print_help();
        exit_status = 1;
    }

    return exit_status;
}
