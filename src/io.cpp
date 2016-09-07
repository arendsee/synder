#include "io.h"

void check_args(size_t line_no, size_t nargs, size_t correct_nargs);

Synmap *load_Synmap(FILE * synfile, int swap, long k, char trans, bool validate)
{
  if(synfile == NULL){
    fprintf(stderr, "NULL synteny input file (%s:%d in %s)", __FILE__, __LINE__, __func__);
    exit(EXIT_FAILURE);
  }

  Synmap *syn = init_Synmap();

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
        SG(syn, loc) = init_Genome(seqid, ncontigs);
      } else if (line[0] == '@') {
        break;
      } else if (line[0] == '$') {
        status = sscanf(line, "$ %zu %zu %s %zu %c\n", &conid, &nblocks, seqid, &contig_length, &dummy);
        check_args(line_no, status, 4);
        SGC(syn, loc, conid) = init_Contig(seqid, nblocks, contig_length, SG(syn, loc));
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
  double score;
  char strand;

  Block *qblk, *tblk;
  Contig *tcon, *qcon;

  while ((read = getline(&line, &len, synfile)) != EOF) {
    line_no++;
    if (line[0] != '$')
      continue;
    unloaded_blocks -= 2;
    status = sscanf(line, "$ %zu %zu %zu %zu %zu %zu %zu %zu %lf %c %c\n",
                    &qcon_id, &qblk_id, &qstart, &qstop,
                    &tcon_id, &tblk_id, &tstart, &tstop, &score, &strand, &dummy);
    check_args(line_no, status, 10);

    qcon = SGC(syn, query, qcon_id);
    tcon = SGC(syn, target, tcon_id);

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
    if (qcon_id >= SG(syn, query)->size
        || tcon_id >= SG(syn, target)->size) {
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

    switch(trans){
        case 'l':
            score = -1 * log(score);
            break;
        case 'd':
            score = score * MIN((tstop - tstart + 1), (qstop - qstart + 1));
            break;
        case 'p':
            score = score * MIN((tstop - tstart + 1), (qstop - qstart + 1)) / 100.0;
            break;
        case 'i':
            // no transformation
            break;
        default:
            fprintf(stderr, "Unexpected transformation '%c'\n", trans);
            exit(EXIT_FAILURE);
            break;
    }

    set_Block(qblk, qstart, qstop, score, '+',    qcon, tblk, line_no);
    set_Block(tblk, tstart, tstop, score, strand, tcon, qblk, line_no);

  }
  free(line);

  /* each block specified has been loaded */
  if (unloaded_blocks != 0) {
    fprintf(stderr, "Unequal number of specified and actual blocks\n");
    exit(EXIT_FAILURE);
  }

  // The following must be run in order
  link_block_corners(syn);
  set_contig_corners(syn);
  merge_all_doubly_overlapping_blocks(syn);
  set_overlap_group(syn);
  link_adjacent_blocks(syn);
  link_contiguous_blocks(syn, k);

  if(validate){
      validate_synmap(syn);
  }

  return (syn);
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
