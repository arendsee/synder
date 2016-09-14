#include "global.h"
#include "arguments.h"
#include "synmap.h"
#include "map.h"

int Offsets::in_start;
int Offsets::in_stop;
int Offsets::out_start;
int Offsets::out_stop;

int main(int argc, char *argv[])
{

    Synmap * syn = NULL;

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

    if (args.synfile)
    {
        syn = new Synmap(args);
        if (args.debug)
        {
            args.print();
            syn->print(true);
            exit(EXIT_SUCCESS);
        }
        if (args.dump_blks)
        {
            syn->dump_blocks();
            exit(EXIT_SUCCESS);
        }
    }

    // No arguments passed
    if (syn == NULL)
    {
        printf("NULL synteny file. Nothing to do ...\n");
        args.print_help();
        exit(EXIT_FAILURE);
    }

    // If no file given by -i, use STDIN, then parse -c options
    if (args.synfile)
    {
        if (args.intfile == NULL)
        {
            args.intfile = stdin;
        }

        if (args.cmd == "")
        {
            fprintf(stderr, "Please, give me a command\n");
        }

        if (args.hitfile)
        {
            fprintf(stderr, "Filter function currently unavailable\n");
        }
        else if (args.cmd == "count")
        {
            syn->count(args.intfile);
        }
        else if (args.cmd == "map")
        {
            syn->map(args.intfile);
        }
        else if (args.cmd == "search")
        {
            find_search_intervals(syn, args.intfile);
        }
        else
        {
            fprintf(stderr, "Command '%s' not recognized\n", args.cmd.c_str());
            args.print_help();
        }
    }

    delete syn;

    return (EXIT_SUCCESS);
}
