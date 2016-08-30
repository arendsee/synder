#ifndef __UI_U__
#define __UI_U__

#include "global.h"

/** Parsed command line arguments */
typedef struct Arguments {
  FILE *synfile;
  FILE *intfile;
  FILE *hitfile;
  char *db_filename;
  char *cmd;
  char **pos;
  char offsets[5];
  long k;
  char trans;
  bool dump_blks;
  bool swap;
  bool debug;
} Arguments;

void close_Arguments(Arguments arg);

Arguments create_Arguments();

void print_version();

void print_help();

void print_args(Arguments args);

Arguments parse_command(int argc, char *argv[]);

#endif
