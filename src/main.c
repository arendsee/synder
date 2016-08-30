#include "global.h"
#include "ui.h"
#include "io.h"
#include "synmap.h"
#include "analysis.h"
#include "map.h"
#include "lev.h"

int main(int argc, char *argv[])
{

  Synmap *syn = NULL;

  // ------------------------------------------------------------------------
  // Prep input 
  // ------------------------------------------------------------------------

  Arguments args = parse_command(argc, argv);

  global_in_start  = (args.offsets[0] == '1');
  global_in_stop   = (args.offsets[1] == '1');
  global_out_start = (args.offsets[2] == '1');
  global_out_stop  = (args.offsets[3] == '1');

  // ------------------------------------------------------------------------
  // Do stuff 
  // ------------------------------------------------------------------------

  // Try to build database. Exit afterwards.
  if (args.db_filename) {
    if (!(args.pos[0] && args.pos[1] && args.pos[2])){
      print_help();
    }

    // TODO: Find a more elegant argument parsing approach
    char db_args[1024];
    sprintf(
        db_args,
        "%s %s -a %s -b %s -i %s -d %s %s %s %s %s 2> /dev/null",
        global_in_start == 1 ? " -x " : "",
        global_in_stop  == 1 ? " -y " : "",
        args.pos[0],
        args.pos[1],
        args.db_filename,
        args.pos[2],
        args.pos[3] != NULL ? "-t" : "",
        args.pos[3] != NULL ? args.pos[3] : "",
        args.pos[4] != NULL ? "-q" : "",
        args.pos[4] != NULL ? args.pos[4] : ""
    );

    bool fail = false;
    char cmd[1024];
    // Try to find the database script
    // Search PATH, working directory, and util folder
    sprintf(cmd, "make-synder-db.sh %s", db_args);
    if(system(cmd) != 0){
        sprintf(cmd, "./make-synder-db.sh %s", db_args);
        if(system(cmd) != 0){
            sprintf(cmd, "util/make-synder-db.sh %s", db_args);
            if(system(cmd) != 0){
                fail = true;
            }
        }
    }
    if(fail){
        fprintf(stderr, "ERROR: Failed to make synder database\n");
        fprintf(stderr, "%s\n", cmd);
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
    syn = load_Synmap(args.synfile, args.swap, args.k, args.trans);
    if(syn != NULL && args.debug){
        print_args(args);
        print_Synmap(syn, true); 
        exit(EXIT_SUCCESS);
    }
  }
  // No arguments passed
  if (syn == NULL) {
    printf("NULL synteny file. Nothing to do ...\n");
    print_help();
  }

  // // Fileter
  // if (args.hitfile) {
  //   if (args.cmd != NULL && strcmp(args.cmd, "filter") == 0) {
  //     size_t width = 5000;
  //     analysis_filter(syn, args.hitfile, single_advocate, &width);
  //   }
  // }
    if (args.hitfile) { printf("Filter function currently unavailable\n"); }

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
      find_search_intervals(syn, args.intfile, false);
    } else if (strcmp(args.cmd, "searchblock") == 0) {
      find_search_intervals(syn, args.intfile, true);
    } else {
      printf("Command '%s' not recognized\n", args.cmd);
      print_help();
    }
  }
  // ------------------------------------------------------------------------
  // Clean up
  // ------------------------------------------------------------------------

  // if (syn != NULL)
  //   free_Synmap(syn);
  close_Arguments(args);

  return (EXIT_SUCCESS);
}
