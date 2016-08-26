#define _GNU_SOURCE
#include <getopt.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


#include "ui.h"
#include "version.h"

#define MAX_POS 5

/** contains all legal arguments from all subcommands */
Arguments create_Arguments()
{
  Arguments args = {
    .synfile     = NULL,
    .intfile     = NULL,
    .hitfile     = NULL,
    .offsets     = "0000",
    .k           = 0,
    .db_filename = NULL,
    .cmd         = NULL,
    .debug       = false,
    .pos         = (char **) malloc(MAX_POS * sizeof(char *))
  };
  memset(args.pos, 0, MAX_POS * sizeof(char *));
  return args;
}

void close_Arguments(Arguments arg)
{
  if (arg.synfile != NULL)
    fclose(arg.synfile);
  if (arg.intfile != NULL)
    fclose(arg.intfile);
  if (arg.hitfile != NULL)
    fclose(arg.hitfile);
  if (arg.db_filename != NULL)
    free(arg.db_filename);
  if (arg.pos != NULL) {
    for (int i = 0; i < 3; i++) {
      if (arg.pos[i] != NULL)
        free(arg.pos[i]);
    }
    free(arg.pos);
  }
  if (arg.cmd != NULL)
    free(arg.cmd);
}

void print_args(Arguments args)
{
  printf(
      "arguments: k=%ld offsets=%s cmd=%s\n",
      args.k,
      args.offsets,
      args.cmd
  );
}

void check_file(FILE * fp, char *name)
{
  if (fp == NULL) {
    fprintf(stderr, "ERROR: Failed to open '%s'\n", name);
    exit(EXIT_FAILURE);
  }
}

void print_version()
{
    printf("%s\n", VERSION);    
    exit(EXIT_SUCCESS);
}

void print_help()
{
  printf("USAGE\n"
         "  synder -d <SYNTENY_FILE> <QUERY> <TARGET> <DB_DIR> <TARGET_GF> <QUERY_GF>\n"
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
         "COMMANDS\n"
         "\tmap\n"
         "\t\t print target intervals overlapping each query interval\n"
         "\tcount\n"
         "\t\t like map but prints only the number that overlap\n"
         "\tfilter\n"
         "\t\t print query-to-target links consistent with the synteny map\n"
         "\tsearch\n"
         "\t\t predict target search spaces for each query interval\n"
         "\tsearchblock\n"
         "\t\t as per search, but also prints blocks\n"
         "\tconvert\n"
         "\t\t convert names in provided gff file to match names in synteny db\n"
         "EXAMPLES\n"
         "  synder -b 1100 -d a-b.syn a b db\n"
         "  synder -b 1111 -i a.gff   -s db/a_b.txt -c map\n"
         "  synder -b 1111 -i a.gff   -s db/a_b.txt -c count\n"
         "  synder -b 1111 -i a.gff   -s db/a_b.txt -c search\n"
         "  synder -b 1111 -i a.gff   -s db/a_b.txt -c convert -t\n"
         "  synder -b 1111 -f hits.syn -s db/a_b.txt -c filter\n"
  );
  exit(EXIT_SUCCESS);
}

Arguments parse_command(int argc, char *argv[])
{
  if (argc == 1)
    print_help();
  int opt;
  FILE *temp;
  Arguments args = create_Arguments();
  while ((opt = getopt(argc, argv, "hvrDd:s:i:c:f:b:k:")) != -1) {
    switch (opt) {
      case 'h':
        print_help();
        break;
      case 'v':
        print_version();
        break;
      case 'D':
        args.debug = true;
        break;
      case 'd':
        temp = fopen(optarg, "r");
        check_file(temp, optarg);
        fclose(temp);
        args.db_filename = strdup(optarg);
        break;
      case 's':
        args.synfile = fopen(optarg, "r");
        check_file(args.synfile, optarg);
        break;
      case 'i':
        args.intfile = fopen(optarg, "r");
        check_file(args.intfile, optarg);
        break;
      case 'f':
        args.hitfile = fopen(optarg, "r");
        check_file(args.hitfile, optarg);
        break;
      case 'c':
        args.cmd = strdup(optarg);
        break;
      case 'r':
        args.swap = true;
        break;
      case 'b':
        strncpy(args.offsets, optarg, 5);
        break;
      case 'k':
        errno = 0;
        long k = strtol(optarg, NULL, 0);
        // see errno docs
        if((errno == ERANGE && (k == LONG_MAX || k == LONG_MIN)) || (errno != 0 && k == 0))
        {
            perror("In argument -k NUM, NUM must be an integer");
        }
        args.k = k;
        break;
      case '?':
        exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; optind < argc; optind++, i++) {
    if (i >= MAX_POS) {
      printf("There can be only %d positionals\n", MAX_POS);
      exit(EXIT_FAILURE);
    }
    args.pos[i] = strdup(argv[optind]);
  }
  return args;
}
