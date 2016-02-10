#ifndef __UI_H__
#define __UI_U__

#include "global.h"

/** Parsed command line arguments */
typedef struct Arguments {
    FILE * synfile;
    FILE * intfile;
    char * db_filename;
    char * cmd;
    char ** pos;
} Arguments;

void close_Arguments(Arguments arg);

Arguments create_Arguments();

void print_help();

void print_args(Arguments args);

Arguments parse_command(int argc, char * argv[]);

#endif
