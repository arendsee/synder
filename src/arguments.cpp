#include "arguments.h"

void Arguments::set_defaults()
{
    synfile     = NULL;
    intfile     = NULL;
    hitfile     = NULL;
    cmd         = "";
    offsets[0]  = 0;
    offsets[1]  = 0;
    offsets[2]  = 0;
    offsets[3]  = 0;
    k           = 0;
    trans       = 'i';
    dump_blks   = false;
    swap        = false;
    debug       = false;
}

Arguments::Arguments(int argc, char *argv[])
{
    set_defaults();

    if (argc == 1)
    {
        print_help();
        exit(EXIT_FAILURE);
    }

    int opt;

    while ((opt = getopt(argc, argv, "hvrDBd:s:i:c:f:b:k:x:")) != -1)
    {
        switch (opt)
        {
        case 'h':
            print_help();
            exit(EXIT_SUCCESS);
            break;
        case 'v':
            print_version();
            exit(EXIT_SUCCESS);
            break;
        case 'D':
            debug = true;
            break;
        case 'B':
            dump_blks = true;
            break;
        case 'x':
            trans = optarg[0];
            if (! (
                        trans == 'i' ||
                        trans == 'd' ||
                        trans == 'p' ||
                        trans == 'l'))
            {
                fprintf(stderr, "-x only takes arguments 'i', 'n', 'l' and 'm'\n");
                exit(EXIT_FAILURE);
            }
            break;
        case 's':
            synfile = fopen(optarg, "r");
            Arguments::check_file(synfile, optarg);
            break;
        case 'i':
            intfile = fopen(optarg, "r");
            Arguments::check_file(intfile, optarg);
            break;
        case 'f':
            hitfile = fopen(optarg, "r");
            Arguments::check_file(hitfile, optarg);
            break;
        case 'c':
            cmd = optarg;
            break;
        case 'r':
            swap = true;
            break;
        case 'b':
            Arguments::set_offsets(optarg);
            break;
        case 'k':
            errno = 0;
            k = strtol(optarg, NULL, 0);
            if ((errno == ERANGE && (k == LONG_MAX || k == LONG_MIN)) || (errno != 0 && k == 0))
            {
                perror("In argument -k NUM, NUM must be an integer");
            }
            break;
        case '?':
            fprintf(stderr, "Argument '%c' not recognized\n", opt);
            exit(EXIT_FAILURE);
        }
    }

    for (; optind < argc; optind++)
    {
        std::string arg = argv[optind];
        pos.push_back(arg);
    }
}


Arguments::~Arguments()
{
    if (synfile != NULL)
        fclose(synfile);
    if (intfile != NULL)
        fclose(intfile);
    if (hitfile != NULL)
        fclose(hitfile);
}


void Arguments::set_offsets(char * optarg)
{
    for (int i = 0; i < 4; i++)
    {
        switch (optarg[i]) {
            case '0':
                offsets[i] = 0;
                break;
            case '1':
                offsets[i] = 1;
                break;
            default:
                fprintf(stderr, "Offsets must be a string of 4 characters, e.g. '0011'\n");
                exit(EXIT_FAILURE);
        }
    }
}


void Arguments::print()
{
    fprintf(
        stderr,
        "arguments: k=%ld offsets=%d%d%d%d cmd=%s\n",
        k,
        offsets[0],
        offsets[1],
        offsets[2],
        offsets[3],
        cmd.c_str()
    );
}


void Arguments::check_file(FILE * fp, char *name)
{
    if (fp == NULL)
    {
        fprintf(stderr, "ERROR: Failed to open '%s'\n", name);
        exit(EXIT_FAILURE);
    }
}


void Arguments::print_version()
{
    printf("%s\n", VERSION);
}


void Arguments::print_help()
{
    printf("USAGE\n"
           "  synder [-t -i <GFF_FILE>] -s <SYNTENY_DB> -c <COMMAND>\n"
           "  synder -f <GFF_FILE> -s <SYNTENY_DB> -c <COMMAND>\n"
           " \n"
           "\t-d \t Convert synteny file into synteny database.\n"
           "\t-f \t Filter using provided gff file.\n"
           "\t-i \t GFF search interval file, if not provided, uses stdin\n"
           "\t-t \t Switch target and query in synteny database\n"
           "\t-c \t Synder command to run. (See Below)\n"
           "\t-b \t Start and stop offsets for input and output (e.g. 1100)\n"
           "\t-k \t Number of interrupting intervals allowed before breaking contiguous set (default=0)\n"
           "\t-D \t Print debug info\n"
           "\t-B \t Dump synteny map links with contiguous set ids\n"
           "\t-x \t Transform score (Synder requires additive scores):\n"
           "\t   \t -'i' := S            (default, no transformation)\n"
           "\t   \t -'d' := L * S        (score densities)\n"
           "\t   \t -'p' := L * S / 100  (percent identity)\n"
           "\t   \t -'l' := -log(S)      (e-values or p-values)\n"
           "\t   \t   Where S is input score and L interval length\n"
           "COMMANDS\n"
           "\tsearch\n"
           "\t\t predict target search spaces for each query interval\n"
           "\tmap\n"
           "\t\t print target intervals overlapping each query interval\n"
           "\tcount\n"
           "\t\t like map but prints only the number that overlap\n"
           "EXAMPLES\n"
           "  synder -b 1100 -i a.gff -s a-b.map -c map\n"
           "  synder -b 1111 -i a.gff -s a-b.map -c count\n"
           "  synder -b 1111 -i a.gff -s a-b.map -c search\n"
          );
}
