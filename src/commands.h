#ifndef __COMMANDS_H__
#define __COMMANDS_H__

#include "optionparser.h"
#include "synmap.h"

#include <vector>
#include <cstring>

bool subcommand_dump   ( int argc, char* argv[] );
bool subcommand_filter ( int argc, char* argv[] );
bool subcommand_map    ( int argc, char* argv[] );
bool subcommand_count  ( int argc, char* argv[] );
bool subcommand_search ( int argc, char* argv[] );

FILE* check_open(const char* name, const char* adjective);
void set_arguments(Arguments& arg, std::vector<option::Option> options);

// Global list of options
enum optionIndex {
    INPUT_GFF,
    INPUT_MAP,
    SYNMAP,
    TCLFILE,
    QCLFILE,
    REVERSE,
    BASE_OFFSETS,
    K,
    TRANSFORM,
    UNKNOWN,
    VERSION,
    HELP // help has to go last
};

#define SHIFT argc -= 1; argv += 1;

#define BOILERPLATE                                                   \
    if (argc == 0) {                                                  \
        option::printUsage(std::cout, usage);                         \
        exit(EXIT_FAILURE);                                           \
    }                                                                 \
                                                                      \
    option::Stats stats(usage, argc, argv);                           \
    std::vector<option::Option> options(stats.options_max);           \
    std::vector<option::Option> buffer(stats.buffer_max);             \
    option::Parser parse(usage, argc, argv, &options[0], &buffer[0]); \
                                                                      \
    if (parse.error()){                                               \
        exit(EXIT_FAILURE);                                           \
    }                                                                 \
                                                                      \
    if (options[HELP]) {                                              \
        option::printUsage(std::cout, usage);                         \
        exit(EXIT_SUCCESS);                                           \
    }                                                                 \

#define SUBCOMMAND_BOILERPLATE                                        \
    BOILERPLATE                                                       \
    Arguments args;                                                   \
    set_arguments(args, options);                                     \
    Synmap syn(args);                                                 \

#endif
