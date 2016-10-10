#include "commands.h"

// Adapted from Lean Mean C++ Parser example scripts
struct Arg: public option::Arg {
    static void printError(const char* msg1, const option::Option& opt, const char* msg2)
    {
        fprintf(stderr, "%s", msg1);
        fwrite(opt.name, opt.namelen, 1, stderr);
        fprintf(stderr, "%s", msg2);
    }

    static option::ArgStatus Required(const option::Option& option, bool msg)
    {
        if (option.arg != 0)
            return option::ARG_OK;

        if (msg) printError("Option '", option, "' requires an argument\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus Filename(const option::Option& option, bool msg)
    {
        if (option.arg != 0) {
            FILE* fp = fopen(option.arg, "r");
            if (fp != nullptr) {
                return option::ARG_OK;
            }
            if (msg) printError("Option '", option, "' file failed to open\n");
        } else {
            if (msg) printError("Option '", option, "' requires an argument\n");
        }
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus Forbidden(const option::Option& option, bool msg)
    {
        if (msg) printError("Option '", option, "' illegal within this subcommand\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus Filename_or_STDIN(const option::Option& option, bool msg)
    {
        if (option.arg != 0) {
            FILE* fp = fopen(option.arg, "r");
            if (fp != nullptr || strcmp(option.arg, "-") == 0) {
                return option::ARG_OK;
            }
            if (msg) printError("Option '", option, "' file failed to open\n");
        } else {
            return option::ARG_OK;
        }
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus Numeric(const option::Option& option, bool msg)
    {
        char* endptr = 0;
        if (option.arg != 0 && strtol(option.arg, &endptr, 10)) {};
        if (endptr != option.arg && *endptr == 0)
            return option::ARG_OK;

        if (msg) printError("Option '", option, "' requires a numeric argument\n");
        return option::ARG_ILLEGAL;
    }
};

namespace descriptors
{
option::Descriptor help = {
    HELP, 0, "h", "help", option::Arg::None,
    "  -h, --help \tPrint usage and exit."
};

option::Descriptor input_gff = {
    INPUT_GFF, 0, "i", "input", Arg::Filename_or_STDIN,
    "  -i, --input \tGFF file of input intervals (default: STDIN)"
};

option::Descriptor input_map = {
    INPUT_MAP, 0, "i", "input", Arg::Filename_or_STDIN,
    "  -i, --input \tA TAB-delimited map file (default: STDIN)\n"
    "              \tThe first 6 fields must be:\n"
    "              \t  1. qseqid - query contig id (e.g. Chr1)\n"
    "              \t  2. qstart - query interval start\n"
    "              \t  3. qstop  - query interval stop\n"
    "              \t  4. sseqid - target contig id\n"
    "              \t  5. sstart - target interval start\n"
    "              \t  6. sstop  - target interval stop\n"
    "              \t"
    "Any number of additional columns may follow, they will be printed without alteration"
};
option::Descriptor synmap = {
    SYNMAP, 0, "s", "synmap", Arg::Filename,
    "  -s, --synmap \tThe primary synteny map (e.g. output of Satsuma)\n"
    "               \tTAB-delimited, no header, with fields\n"
    "               \t  1. qseqid - query contig id (e.g. Chr1)\n"
    "               \t  2. qstart - query interval start\n"
    "               \t  3. qstop  - query interval stop\n"
    "               \t  4. sseqid - target contig id\n"
    "               \t  5. sstart - target interval start\n"
    "               \t  6. sstop  - target interval stop\n"
    "               \t  7. score  - score of the syntenic match*\n"
    "               \t  8. strand - relative orientation\n"
    "               \t  * score can be any numeric value, it will be\n"
    "               \t    transformed as specified by the -x option"
};

option::Descriptor tclfile = {
    TCLFILE, 0, "t", "tcl", Arg::Filename,
    "  -t, --tcl \tTarget chromosome lengths file\n"
    "            \tTAB-delimited with columns: <name>, <length>"
};

option::Descriptor qclfile = {
    QCLFILE, 0, "q", "qcl", Arg::Filename,
    "  -q, --qcl \tQuery chromosome lengths file\n"
    "            \tTAB-delimited with columns: <name>, <length>"
};

option::Descriptor reverse = {
    REVERSE, 0, "r", "reverse", option::Arg::None,
    "  -r, --reverse \tSwitch target and query in synteny map"
};

option::Descriptor base_offsets = {
    BASE_OFFSETS, 0, "b", "base-offsets", Arg::Required,
    "  -b, --base-offsets \tStart and stop offsets for input and output (e.g. 1100)"
};
option::Descriptor forbid_base_offsets = {
    BASE_OFFSETS, 0, "", "", Arg::Forbidden, ""
};

option::Descriptor k = {
    K, 0, "k", "--interruption-fuzz", Arg::Numeric,
    "  -k, --interruption-fuzz \tNumber of interrupting intervals allowed before breaking contiguous set (default=0)"
};

option::Descriptor transform = {
    TRANSFORM, 0, "x", "--transform", Arg::Required,
    "  -x, --transform \tTransform score (Synder requires additive scores):\n"
    "                 \t -'i' := S            (default, no transformation)\n"
    "                 \t -'d' := L * S        (transform from score densities)\n"
    "                 \t -'p' := L * S / 100  (transform from percent identity)\n"
    "                 \t -'l' := -log(S)      (transform from e-values or p-values)\n"
    "                 \t   Where S is input score and L interval length"
};
}

bool subcommand_dump(int argc, char* argv[])
{
    const option::Descriptor usage[] = {
        {
            UNKNOWN, 0, "", "", option::Arg::None,
            "synder dump - print all blocks with contiguous set ids\n"
            "USAGE\n"
            "  synder dump -s SYNTENY_MAP\n"
            "DESCRIPTION\n"
            "  Dump processed blocks to stdout and exit. The output contains\n"
            "  all the fields in the synteny map, with preserved order, and\n"
            "  adds a column of contiguous set ids. Also, doubly overlapping\n"
            "  bocks are also merged.\n"
            "ARGUMENTS"
        },
        descriptors::synmap,
        descriptors::reverse,
        descriptors::transform,
        descriptors::base_offsets,
        descriptors::help,
        {0,0,0,0,0,0}
    };

    try {
        SUBCOMMAND_BOILERPLATE
        syn.dump_blocks();
        return true;
    } catch (...) {
        return false;
    }
}

bool subcommand_filter(int argc, char* argv[])
{
    const option::Descriptor usage[] = {
        {
            UNKNOWN, 0, "", "", option::Arg::None,
            "synder filter - remove links that disagree with the synteny map\n"
            "USAGE\n"
            "  synder filter -i INPUT.gff -s SYNMAP.tab [OPTIONS]\n"
            "DESCRIPTION\n"
            "  Remove input links that disagree with the synteny map Given a map\n"
            "  between the target and query (e.g. BLAST output) and a synteny map,\n"
            "  remove all entries in the map that do not overlap search intervals\n"
            "  in the query.\n"
            "ARGUMENTS"
        },
        descriptors::synmap,
        descriptors::input_map,
        descriptors::tclfile,
        descriptors::qclfile,
        descriptors::k,
        descriptors::reverse,
        descriptors::help,
        {0,0,0,0,0,0}
    };

    try {
        SUBCOMMAND_BOILERPLATE
        syn.filter(args.intfile);;
        return true;
    } catch (...) {
        return false;
    }
}

bool subcommand_map(int argc, char* argv[])
{
    const option::Descriptor usage[] = {
        {
            UNKNOWN, 0, "", "", option::Arg::None,
            "synder map - trace intervals across genomes\n"
            "USAGE\n"
            "  synder map -i INPUT.gff -s SYNMAP.tab [OPTIONS]\n"
            "DESCRIPTION\n"
            "  Given a set of query intervals, find all synmap queries that\n"
            "  overlap and map to homologous target intervals.\n"
            "ARGUMENTS"
        },
        descriptors::synmap,
        descriptors::input_gff,
        descriptors::reverse,
        descriptors::base_offsets,
        descriptors::help,
        {0,0,0,0,0,0}
    };

    try {
        SUBCOMMAND_BOILERPLATE
        syn.process_gff(args.intfile, C_MAP);;
        return true;
    } catch (...) {
        return false;
    }
}

bool subcommand_count(int argc, char* argv[])
{
    const option::Descriptor usage[] = {
        {
            UNKNOWN, 0, "", "", option::Arg::None,
            "synder count - count overlaps\n"
            "USAGE\n"
            "  synder count -i INPUT.gff -s SYNMAP.tab [OPTIONS]\n"
            "DESCRIPTION\n"
            "  Given a set of query intervals, count all synmap queries that\n"
            "  overlap and map to homologous target intervals.\n"
            "ARGUMENTS"
        },
        descriptors::synmap,
        descriptors::input_gff,
        descriptors::reverse,
        descriptors::base_offsets,
        descriptors::help,
        {0,0,0,0,0,0}
    };

    try {
        SUBCOMMAND_BOILERPLATE
        syn.process_gff(args.intfile, C_COUNT);;
        return true;
    } catch (...) {
        return false;
    }
}

bool subcommand_search(int argc, char* argv[])
{
    const option::Descriptor usage[] = {
        {
            UNKNOWN, 0, "", "", option::Arg::None,
            "synder search - predict search intervals\n"
            "USAGE\n"
            "  synder search -i INPUT.gff -s SYNMAP.tab [OPTIONS]\n"
            "DESCRIPTION\n"
            "  Given an input GFF and a synteny, find search intervals. See README for\n"
            "  algorithmic details and a discussion of scoring and argument -k.\n"
            "ARGUMENTS"
        },
        descriptors::synmap,
        descriptors::input_gff,
        descriptors::tclfile,
        descriptors::qclfile,
        descriptors::k,
        descriptors::reverse,
        descriptors::base_offsets,
        descriptors::transform,
        descriptors::help,
        {
            UNKNOWN, 0, "", "", option::Arg::None,
            "EXAMPLE\n"
            "  # Using files in the sample-inputs/ folder\n"
            "  synder search -i a.gff -s a-b.tab -b 0011 -x p\n"
            "  # Search the opposite direction\n"
            "  synder search -r -i b.gff -s a-b.tab\n"
        },
        {0,0,0,0,0,0}
    };

    try {
        SUBCOMMAND_BOILERPLATE
        syn.process_gff(args.intfile, C_SEARCH);;
        return true;
    } catch (...) {
        return false;
    }
}

void set_arguments(Arguments& arg, std::vector<option::Option> options)
{
    if(options[INPUT_GFF].desc != 0) {
        arg.intfile = check_open(options[INPUT_GFF].arg, "GFF file");
    }

    if(options[INPUT_MAP].desc != 0) {
        arg.intfile = check_open(options[INPUT_MAP].arg, "input map");
    }

    if(options[SYNMAP].desc != 0) {
        arg.synfile = check_open(options[SYNMAP].arg, "synteny map");
    }

    if(options[TCLFILE].desc != 0) {
        arg.tclfile = check_open(options[TCLFILE].arg, "target contig length file");
    }

    if(options[QCLFILE].desc != 0) {
        arg.qclfile = check_open(options[QCLFILE].arg, "query contig length file");
    }

    if(options[REVERSE].desc != 0) {
        arg.swap = true;
    }

    if(options[BASE_OFFSETS].desc != 0) {
        std::string a(options[BASE_OFFSETS].arg);
        if(a.length() == 4) {
            Offsets::in_start  = a[0] == '0' ? 0 : 1;
            Offsets::in_stop   = a[1] == '0' ? 0 : 1;
            Offsets::out_start = a[2] == '0' ? 0 : 1;
            Offsets::out_stop  = a[3] == '0' ? 0 : 1;
        } else {
            fprintf(stderr, "Offsets must be a string of length 4\n");
            exit(EXIT_FAILURE);
        }
    }

    if(options[K].desc != 0) {
        long k = strtol(options[K].arg, nullptr, 0);
        if ((errno == ERANGE && (k == LONG_MAX || k == LONG_MIN)) || (errno != 0 && k == 0)) {
            perror("In argument -k NUM, NUM must be an integer");
        }
        arg.k = k;
    }

    if(options[TRANSFORM].desc != 0) {
        char trans = options[TRANSFORM].arg[0];
        if (! (trans == 'i' || trans == 'd' || trans == 'p' || trans == 'l')) {
            fprintf(stderr, "-x only takes arguments 'i', 'd', 'p' and 'l'\n");
            exit(EXIT_FAILURE);
        }
        arg.trans = trans;
    }

    // If no file given by -i, use STDIN
    if (arg.intfile == nullptr) {
        arg.intfile = stdin;
    }
}

FILE* check_open(const char* name, const char* adjective)
{
    FILE* fp = fopen(name, "r");
    if (fp == nullptr) {
        fprintf(stderr, "ERROR: Failed to open %s '%s'\n", adjective, name);
        exit(EXIT_FAILURE);
    }
    return fp;
}
