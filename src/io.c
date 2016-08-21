#include <assert.h>

#include "io.h"

void check_args(size_t line_no, size_t nargs, size_t correct_nargs);

Synmap *load_Synmap(FILE * synfile, int swap)
{
  assert(synfile != NULL);

  Synmap *synmap = init_Synmap();

  int query = swap;
  int target = !swap;
  size_t line_no = 0;
  size_t unloaded_blocks = 0;
  int status;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  char seqid[128];
  size_t ncontigs;
  size_t nblocks;
  size_t conid;
  size_t contig_length;

  // NOTE: The synteny database files are always 0-based. So subtracting
  // global_in_base from the intervals is not necessary.

  // A dummy variable searched for at the end of each sscanf format sequence.
  // It will be matched only if there are too many arguments.
  char dummy;
  for (int i = 0; i < 2; i++) {
    while ((read = getline(&line, &len, synfile)) != EOF) {
      line_no++;
      size_t loc = i == query ? 0 : 1;
      if (line[0] == '>') {
        status = sscanf(line, "> %s %zu %c", seqid, &ncontigs, &dummy);
        check_args(line_no, status, 2);
        SG(synmap, loc) = init_Genome(seqid, ncontigs);
      } else if (line[0] == '@') {
        break;
      } else if (line[0] == '$') {
        status = sscanf(line, "$ %zu %zu %s %zu %c\n", &conid, &nblocks, seqid, &contig_length, &dummy);
        check_args(line_no, status, 4);
        SGC(synmap, loc, conid) = init_Contig(seqid, nblocks, contig_length);
        unloaded_blocks += nblocks;
      } else {
        fprintf(stderr, "Incorrect file format, line %zu\n", line_no);
        fprintf(stderr, "Offending line:\n%s", line);
        exit(EXIT_FAILURE);
      }
    }
  }


  line_no = 0;
  size_t qcon_id, qblk_id, tcon_id, tblk_id;
  long qstart, qstop, tstart, tstop;
  float score;
  char strand;

  Block *qblk, *tblk;
  Contig *tcon, *qcon;

  while ((read = getline(&line, &len, synfile)) != EOF) {
    line_no++;
    if (line[0] != '$')
      continue;
    unloaded_blocks -= 2;
    status = sscanf(line, "$ %zu %zu %zu %zu %zu %zu %zu %zu %f %c %c\n",
                    &qcon_id, &qblk_id, &qstart, &qstop,
                    &tcon_id, &tblk_id, &tstart, &tstop, &score, &strand, &dummy);
    check_args(line_no, status, 10);

    qcon = SGC(synmap, query, qcon_id);
    tcon = SGC(synmap, target, tcon_id);

    if (qstart > qstop || tstart > tstop) {
      fprintf(stderr, "start must be less than stop on line %zu\n", line_no);
      fprintf(stderr, "offending line:\n%s\n", line);
      exit(EXIT_FAILURE);
    }
    if (tcon->length <= tstop) {
      fprintf(stderr, "stop must be less than contig length on line %zu\n", line_no);
      fprintf(stderr, "conid=%zu blkid=%zu pos=(%zu, %zu) conlen=%zu\n",
        tcon_id, tblk_id, tstart, tstop, tcon->length); 
      fprintf(stderr, "offending line:\n%s\n", line);
      exit(EXIT_FAILURE);
    }
    // don't exceed the specified number of Contig in Genome
    if (qcon_id >= SG(synmap, query)->size
        || tcon_id >= SG(synmap, target)->size) {
      fprintf(stderr, "too few contigs specified\n");
      exit(EXIT_FAILURE);
    }
    // don't exceed the specified size of the Contig block arrays
    if (qblk_id >= qcon->size ||
        tblk_id >= tcon->size) {
      fprintf(stderr, "too few blocks specified\n");
      exit(EXIT_FAILURE);
    }

    qblk = &qcon->block[qblk_id];
    tblk = &tcon->block[tblk_id];

    set_Block(qblk, qstart, qstop, score, '+',    qcon, tblk);
    set_Block(tblk, tstart, tstop, score, strand, tcon, qblk);

  }
  free(line);

  /* each block specified has been loaded */
  if (unloaded_blocks != 0) {
    fprintf(stderr, "Unequal number of specified and actual blocks\n");
    exit(EXIT_FAILURE);
  }

  link_four_corners(synmap);

  set_head_and_tail(synmap);

  set_overlap_group(synmap);

  link_adjacent_blocks(synmap);

  link_contiguous_blocks(synmap);

  validate_synmap(synmap);

  return (synmap);
}

void check_args(size_t line_no, size_t nargs, size_t correct_nargs)
{
  if (nargs < correct_nargs) {
    fprintf(stderr, "Either too few fields on input line %zu,", line_no);
    fprintf(stderr, "or they are of the wrong type\n");
    fprintf(stderr, "Found %zu arguments, expected %zu.\n", nargs, correct_nargs);
    exit(EXIT_FAILURE);
  } else if (nargs > correct_nargs) {
    fprintf(stderr, "Too many fields on input line %zu\n", line_no);
    exit(EXIT_FAILURE);
  }
}
