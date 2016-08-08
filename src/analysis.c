#include "analysis.h"

void analysis_count(Synmap * syn, FILE * intfile)
{
  char seqname[128];
  uint count;
  int chrid, start, stop;
  while ((fscanf(intfile,
                 "%d %*s %*s %d %d %*s %*c %*s %s\n",
                 &chrid, &start, &stop, seqname)) != EOF) {
    // contig.c::count_overlaps ->
    //   itree/search.c::count_interval_overlaps ->
    //   itree/search.c::count_interval_overlaps_r
    count = count_overlaps(SGC(syn, 0, chrid), start, stop);
    printf("%s\t%u\n", seqname, count);
  }
}

void analysis_map(Synmap * syn, FILE * intfile)
{
  char seqname[128];
  int chrid, start, stop;
  ResultContig * rc;
  Contig *contigs;
  Contig *tcon;
  Block *qblk, *tblk;
  bool missing;
  while ((fscanf(intfile,
                 "%d %*s %*s %d %d %*s %*c %*s %s\n",
                 &chrid, &start, &stop, seqname)) != EOF) {

    rc = get_region(SGC(syn, 0, chrid), start, stop);
    contigs = rc->contig;
    missing = rc->inbetween || rc->leftmost || rc->rightmost;

    for (int i = 0; i < contigs->size; i++) {
      qblk = contigs->block[i];
      if (qblk != NULL) {
        tcon = qblk->over->parent;
        tblk = qblk->over;
        printf("%s %s %u %u %d\n",
               seqname, tcon->name, tblk->pos[0], tblk->pos[1], missing);
      }
    }

    free(rc);
    free(contigs->name);
    free(contigs->block);
    free(contigs);
  }
}
