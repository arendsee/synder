#include "analysis.h"

void analysis_count(Synmap * syn, FILE * intfile)
{
  char seqname[128];
  size_t count;
  size_t chrid;
  long start, stop;
  while ((fscanf(intfile,
                 "%zu %*s %*s %li %li %*s %*c %*s %s\n",
                 &chrid, &start, &stop, seqname)) != EOF)
  {
    check_in_offset(start, stop);
    start -= global_in_start;
    stop  -= global_in_stop;

    // contig.c::count_overlaps ->
    //   itree/search.c::count_interval_overlaps ->
    //   itree/search.c::count_interval_overlaps_r
    count = count_overlaps(SGC(syn, 0, chrid), start, stop);
    printf("%s\t%zu\n", seqname, count);
  }
}

void analysis_map(Synmap * syn, FILE * intfile)
{
  char seqname[128];
  size_t chrid;
  long start, stop;
  ResultContig * rc;
  Block *qblk, *tblk;
  bool missing;
  while ((fscanf(intfile,
                 "%zu %*s %*s %zu %zu %*s %*c %*s %s\n",
                 &chrid, &start, &stop, seqname)) != EOF)
  {
    check_in_offset(start, stop);
    start -= global_in_start;
    stop  -= global_in_stop;

    rc = get_region(SGC(syn, 0, chrid), start, stop, false);
    missing = rc->inbetween || rc->leftmost || rc->rightmost;

    for (size_t i = 0; i < rc->size; i++) {
      qblk = rc->block[i];
      if (qblk != NULL) {
        tblk = qblk->over;
        printf("%s %s %zu %zu %d\n",
               seqname,
               tblk->parent->name,
               tblk->pos[0] + global_out_start,
               tblk->pos[1] + global_out_stop,
               missing
        );
      }
    }

    free(rc);
  }
}
