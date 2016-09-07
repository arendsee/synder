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
    start -= Offsets::in_start;
    stop  -= Offsets::in_stop;

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
    start -= Offsets::in_start;
    stop  -= Offsets::in_stop;

    rc = get_region(SGC(syn, 0, chrid), start, stop, false);
    missing = rc->inbetween || rc->leftmost || rc->rightmost;

    for (size_t i = 0; i < rc->size; i++) {
      qblk = rc->block[i];
      if (qblk != NULL) {
        tblk = qblk->over;
        printf("%s %s %zu %zu %d\n",
               seqname,
               tblk->parent->name,
               tblk->pos[0] + Offsets::out_start,
               tblk->pos[1] + Offsets::out_stop,
               missing
        );
      }
    }

    free(rc);
  }
}
