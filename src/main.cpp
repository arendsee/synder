#include "global.h"
#include "arguments.h"
#include "synmap.h"

int main(int argc, char *argv[])
{

    int exit_status = 0;

    // ------------------------------------------------------------------------
    // Prep input
    // ------------------------------------------------------------------------

    Arguments args = Arguments(argc, argv);

    Offsets::in_start  = args.offsets[0];
    Offsets::in_start  = args.offsets[1];
    Offsets::out_start = args.offsets[2];
    Offsets::out_start = args.offsets[3];

    // ------------------------------------------------------------------------
    // Do stuff
    // ------------------------------------------------------------------------

    // If no file given by -i, use STDIN
    if (args.intfile == nullptr)
        args.intfile = stdin;

    if (args.synfile) {
        Synmap syn(args);
        switch(args.cmd) {
            case C_DEBUG:
                args.print();
                syn.print(true);
                break;
            case C_DUMP:
                syn.dump_blocks();
                break;
            default:
                exit_status = syn.process_gff(args.intfile, args.cmd);
        }
    }
    else {
        printf("Nothing to do ...\n");
        args.print_help();
        exit_status = 1;
    }

    return exit_status;
}
