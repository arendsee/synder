#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "ui.h"
#include "io.h"
#include "synmap.h"
#include "analysis.h"
#include "contiguous.h"
#include "lev.h"

int main(int argc, char *argv[])
{


  Synmap *syn = NULL;

  // ------------------------------------------------------------------------
  // Prep input 
  // ------------------------------------------------------------------------

  Arguments args = parse_command(argc, argv);

  // ------------------------------------------------------------------------
  // Do stuff 
  // ------------------------------------------------------------------------

  // Try to build database. Exit afterwards.
  if (args.db_filename) {
    if (!(args.pos[0] && args.pos[1] && args.pos[2])){
      print_help();
    }
    char cmd[512];
    sprintf(cmd,
            "make-synder-db.sh -a %s -b %s -i %s -d %s",
            args.pos[0], args.pos[1], args.db_filename, args.pos[2]);
    
    int exit_status = system(cmd);
    if(exit_status != 0){
        fprintf(stderr, "ERROR: Failed to make synder database\n");
        exit(EXIT_FAILURE);
    } else {
        exit(EXIT_SUCCESS);
    }
  }
  // Converts between gff contig naming conventions (field 1)
  if (args.cmd != NULL && strcmp(args.cmd, "convert") == 0 && args.intfile && args.intfile) {
    convert_seqname(args.synfile, args.intfile, args.swap);
    exit(EXIT_SUCCESS);
  }
  // Load synteny db 
  if (args.synfile) {
    syn = load_Synmap(args.synfile, args.swap);
  }
  // No arguments passed
  if (syn == NULL) {
    printf("NULL synteny file. Nothing to do ...\n");
    print_help();
  }
  // Fileter
  if (args.hitfile) {
    if (args.cmd != NULL && strcmp(args.cmd, "filter") == 0) {
      int width = 5000;
      analysis_filter(syn, args.hitfile, single_advocate, &width);
    }
  }
  // If no file given by -i, use STDIN, then parse -c options
  if (args.synfile) {
    if (args.intfile == NULL) {
      args.intfile = stdin;
    }

    if (args.cmd == NULL) {
      printf("Please, give me a command\n");
    } else if (strcmp(args.cmd, "count") == 0) {
      analysis_count(syn, args.intfile);
    } else if (strcmp(args.cmd, "map") == 0) {
      analysis_map(syn, args.intfile);
    } else if (strcmp(args.cmd, "search") == 0) {
      contiguous_query(syn, args.intfile, false);
    } else if (strcmp(args.cmd, "searchblock") == 0) {
      contiguous_query(syn, args.intfile, true);
    } else {
      printf("Command '%s' not recognized\n", args.cmd);
      print_help();
    }
  }
  // ------------------------------------------------------------------------
  // Clean up
  // ------------------------------------------------------------------------

  if (syn != NULL)
    free_Synmap(syn);
  close_Arguments(args);

  return (EXIT_SUCCESS);
}
