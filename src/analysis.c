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
  Block *qblk;
  Block *tblk;
  Block twoblk;
  bool missing;
  while ((fscanf(intfile,
                 "%d %*s %*s %d %d %*s %*c %*s %s\n",
                 &chrid, &start, &stop, seqname)) != EOF) {
    missing = false;
    rc = get_region(SGC(syn, 0, chrid), start, stop);
    contigs = rc->contig;
    free(rc);
    // If the interval is between blocks, the size will ALWAYS be 2,
    // However, one of these may be NULL
    if (contigs->size == 2) {
      twoblk.start = start;
      twoblk.stop = stop;
      if (!((CB(contigs, 0) && block_overlap(CB(contigs, 0), &twoblk)) ||
            (CB(contigs, 1) && block_overlap(CB(contigs, 1), &twoblk))
          )) {
        missing = true;
      }
    }
    for (int i = 0; i < contigs->size; i++) {
      qblk = contigs->block[i];
      if (qblk != NULL) {
        tcon = QT_SGC(syn, qblk);
        tblk = QT_SGCB(syn, qblk);
        printf("%s %s %u %u %d\n",
               seqname, tcon->name, tblk->start, tblk->stop, missing);
      }
    }
    free(contigs->name);
    free(contigs->block);
    free(contigs);
  }
}

void analysis_filter(Synmap * syn, FILE * hitfile,
                     bool(*classifier) (Synmap *, Link *, void *), void *arg)
{
  Link link;
  char *line = NULL;
  size_t len = 0;
  int read = 0;
  bool agrees;
  while ((read = getline(&line, &len, hitfile)) != EOF) {
    sscanf(line, "%lu %u %u %lu %u %u\n",
           &link.qseqid, &link.qstart, &link.qstop,
           &link.tseqid, &link.tstart, &link.tstop);
    agrees = classifier(syn, &link, arg);
    printf("%d\t%s", agrees, line);
  }
  free(line);
}

bool single_advocate(Synmap * syn, Link * query, void *width_ptr)
{
  uint width;
  ResultContig * rc;
  Contig *con, *qcon;
  Block *qblk, *tblk;

  width = *(uint *) width_ptr;

  // get the region that actually overlaps (or flanks, if the query is inbetween hits)
  rc = get_region(SGC(syn, 0, query->qseqid), query->qstart, query->qstop);
  con = rc->contig;
  free(rc);

  // determine the search interval on the query
  Block qgood = {
    .start = (query->qstart > width) ? query->qstart - width : 0,
    .stop = query->qstop + width
  };

  // determine the search interval on the target
  Block tgood = {
    .start = (query->tstart > width) ? query->tstart - width : 0,
    .stop = query->tstop + width
  };


  // Query contig pointer
  qcon = SGC(syn, 0, query->qseqid);

  if (con->block[0]) {
    // look down
    int lo_id = CB_STOPID(con, 0);
    for (; lo_id >= 0 && CB_STOP(qcon, lo_id) > qgood.start; lo_id--) {
      qblk = CB(qcon, lo_id);
      if (qblk->oseqid == query->tseqid) {
        tblk = QT_SGCB(syn, qblk);
        if (block_overlap(&tgood, tblk))
          return true;
      }
    }
  }

  if (con->block[con->size - 1]) {
    // look up 
    int hi_id = CB_STARTID(con, con->size - 1);
    for (; hi_id < qcon->size && CB_START(qcon, hi_id) < qgood.stop; hi_id++) {
      qblk = CB(qcon, hi_id);
      if (qblk->oseqid == query->tseqid) {
        tblk = QT_SGCB(syn, qblk);
        if (block_overlap(&tgood, tblk))
          return true;
      }
    }
  }

    /** \todo find the memory leak  */
  free(con->block);
  free(con->name);
  free(con->by_stop);
  free(con);
  return false;
}
