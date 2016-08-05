#include <assert.h>

#include "io.h"

void check_args(int line_no, int nargs, int correct_nargs);

Synmap *load_Synmap(FILE * synfile, int swap)
{
  assert(synfile != NULL);

  Synmap *synmap = init_Synmap();

  int query = swap;
  int target = !swap;
  int line_no = 0;
  int unloaded_blocks = 0;
  int status;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  char seqid[128];
  uint ncontigs, nblocks, id;

  // A dummy variable searched for at the end of each sscanf format sequence.
  // It will be matched only if there are too many arguments.
  char dummy;
  for (int i = 0; i < 2; i++) {
    while ((read = getline(&line, &len, synfile)) != EOF) {
      line_no++;
      int loc = i == query ? 0 : 1;
      if (line[0] == '>') {
        status = sscanf(line, "> %s %u %c", seqid, &ncontigs, &dummy);
        check_args(line_no, status, 2);
        SG(synmap, loc) = init_Genome(seqid, ncontigs);
      } else if (line[0] == '@') {
        break;
      } else if (line[0] == '$') {
        status = sscanf(line, "$ %u %u %s %c\n", &id, &nblocks, seqid, &dummy);
        check_args(line_no, status, 3);
        SGC(synmap, loc, id) = init_Contig(seqid, nblocks);
        unloaded_blocks += nblocks;
      } else {
        fprintf(stderr, "Incorrect file format, line %d\n", line_no);
        fprintf(stderr, "Offending line:\n%s", line);
        exit(EXIT_FAILURE);
      }
    }
  }


  line_no = 0;
  uint qcon_id, qblk_id, qstart, qstop;
  uint tcon_id, tblk_id, tstart, tstop;
  char strand;

  while ((read = getline(&line, &len, synfile)) != EOF) {
    line_no++;
    if (line[0] != '$')
      continue;
    unloaded_blocks -= 2;
    status = sscanf(line, "$ %u %u %u %u %u %u %u %u %c %c\n",
                    &qcon_id, &qblk_id, &qstart, &qstop,
                    &tcon_id, &tblk_id, &tstart, &tstop, &strand, &dummy);
    check_args(line_no, status, 10);
    if (qstart > qstop || tstart > tstop) {
      fprintf(stderr, "start must be less than stop on line %d\n", line_no);
      exit(EXIT_FAILURE);
    }
    // don't exceed the specified number of Contig in Genome
    if (qcon_id >= SG(synmap, query)->size
        || tcon_id >= SG(synmap, target)->size) {
      fprintf(stderr, "too few contigs specified\n");
      exit(EXIT_FAILURE);
    }
    // don't exceed the specified size of the Contig block arrays
    if (qblk_id >= SGC(synmap, query, qcon_id)->size ||
        tblk_id >= SGC(synmap, target, tcon_id)->size) {
      fprintf(stderr, "too few blocks specified\n");
      exit(EXIT_FAILURE);
    }

    SGCB(synmap, query, qcon_id, qblk_id) =
      init_Block(qstart, qstop, tcon_id, tblk_id, strand);

    SGCB(synmap, target, tcon_id, tblk_id) =
      init_Block(tstart, tstop, qcon_id, qblk_id, strand);
  }
  free(line);

  /* each block specified has been loaded */
  if (unloaded_blocks != 0) {
    fprintf(stderr, "Unequal number of specified and actual blocks\n");
    exit(EXIT_FAILURE);
  }

  sort_all_contigs(synmap);

  // Set linkids in ascending order relative to query
  uint linkid = 0; 
  uint gsize, csize;
  for(int g = 0; g < 2; g++){
    gsize = SG(synmap, 0)->size;
    for(int i = 0; i < gsize; i++){
      csize = SGC(synmap, 0, i)->size;
      for(int j = 0; j < csize; j++){
        SGCB(synmap, 0, i, j)->linkid = linkid;
        linkid++;
      }
    }
  }

  return (synmap);
}

void check_args(int line_no, int nargs, int correct_nargs)
{
  if (nargs < correct_nargs) {
    fprintf(stderr, "Either too few fields on input line %d,", line_no);
    fprintf(stderr, "or they are of the wrong type\n");
    fprintf(stderr, "Found %d arguments, expected %d.\n", nargs, correct_nargs);
    exit(EXIT_FAILURE);
  } else if (nargs > correct_nargs) {
    fprintf(stderr, "Too many fields on input line %d\n", line_no);
    exit(EXIT_FAILURE);
  }
}
