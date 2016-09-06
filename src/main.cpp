#include "global.hpp"
#include "arguments.hpp"
#include "io.hpp"
#include "synmap.hpp"
#include "analysis.hpp"
#include "map.hpp"
#include "lev.hpp"

int Offsets::in_start;
int Offsets::in_stop;
int Offsets::out_start;
int Offsets::out_stop;

int main(int argc, char *argv[])
{

    Synmap *syn = NULL;

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

    // Try to build database. Exit afterwards.
    if (args.db_filename != "")
    {
        if (args.pos.size() < 3)
        {
            fprintf(stderr, "Too few positional arguments, 3-5 expected, %zu found\n", args.pos.size()); 
            args.print_help();
            exit(EXIT_FAILURE);
        }

        // TODO: Find a more elegant argument parsing approach
        char db_args[1024];
        sprintf(
            db_args,
            "%s %s -a %s -b %s -i %s -d %s %s %s %s %s 2> /dev/null",
            Offsets::in_start == 1 ? " -x " : "",
            Offsets::in_stop  == 1 ? " -y " : "",
            args.pos.at(0).c_str(),
            args.pos.at(1).c_str(),
            args.db_filename.c_str(),
            args.pos.at(2).c_str(),
            args.pos.size() > 3 ? "-t" : "",
            args.pos.size() > 3 ? args.pos[3].c_str() : "",
            args.pos.size() > 4 ? "-q" : "",
            args.pos.size() > 4 ? args.pos[4].c_str() : ""
        );

        bool fail = false;
        char cmd[1024];
        // Try to find the database script
        // Search PATH, working directory, and util folder
        sprintf(cmd, "make-synder-db.sh %s", db_args);
        if (system(cmd) != 0)
        {
            sprintf(cmd, "./make-synder-db.sh %s", db_args);
            if (system(cmd) != 0)
            {
                sprintf(cmd, "util/make-synder-db.sh %s", db_args);
                if (system(cmd) != 0)
                {
                    fail = true;
                }
            }
        }
        if (fail)
        {
            fprintf(stderr, "ERROR: Failed to make synder database\n");
            fprintf(stderr, "%s\n", cmd);
            exit(EXIT_FAILURE);
        }
        else
        {
            exit(EXIT_SUCCESS);
        }
    }
    // Converts between gff contig naming conventions (field 1)
    if (args.cmd == "convert" && args.intfile && args.intfile)
    {
        convert_seqname(args.synfile, args.intfile, args.swap);
        exit(EXIT_SUCCESS);
    }
    // Load synteny db
    if (args.synfile)
    {
        syn = load_Synmap(args.synfile, args.swap, args.k, args.trans, args.validate);
        if (syn != NULL && args.debug)
        {
            args.print();
            print_Synmap(syn, true);
            exit(EXIT_SUCCESS);
        }
        if (syn != NULL && args.dump_blks)
        {
            dump_blocks(syn);
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

    // // Fileter
    // if (args.hitfile) {
    //   if (args.cmd != NULL && args.cmd == "filter") {
    //     size_t width = 5000;
    //     analysis_filter(syn, args.hitfile, single_advocate, &width);
    //   }
    // }
    if (args.hitfile)
    {
        printf("Filter function currently unavailable\n");
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
            printf("Please, give me a command\n");
        }
        else if (args.cmd == "count")
        {
            analysis_count(syn, args.intfile);
        }
        else if (args.cmd == "map")
        {
            analysis_map(syn, args.intfile);
        }
        else if (args.cmd == "search")
        {
            find_search_intervals(syn, args.intfile);
        }
        else
        {
            printf("Command '%s' not recognized\n", args.cmd.c_str());
            args.print_help();
        }
    }

    return (EXIT_SUCCESS);
}
