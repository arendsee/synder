#include "analysis.h"

void analysis_count(Synmap * syn, FILE * intfile){
    char seqname[128];
    uint count;
    int chrid, start, stop;
    while ((fscanf(intfile,
                   "%d %*s %*s %d %d %*c %*c %*s %s\n",
                   &chrid, &start, &stop, seqname)) != EOF)
    {
        count = count_overlaps(SGC(syn, 0, chrid), start, stop);
        printf("%s\t%u\n", seqname, count);
    }
}

void analysis_map(Synmap * syn, FILE * intfile){
    char seqname[128];
    int chrid, start, stop;
    Contig * contigs;
    Contig * tcon;
    Block * qblk;
    Block * tblk;
    Block twoblk;
    bool missing;
    while ((fscanf(intfile,
                   "%d %*s %*s %d %d %*c %*c %*s %s\n",
                   &chrid, &start, &stop, seqname)) != EOF)
    {
        missing = false;
        contigs = get_region(SGC(syn, 0, chrid), start, stop);
        // If the interval is between blocks, the size will ALWAYS be 2,
        // However, one of these may be NULL
        if(contigs->size == 2){
            twoblk.start = start;
            twoblk.stop = stop;
            if(!((CB(contigs, 0) && block_overlap(CB(contigs, 0), &twoblk)) || 
                 (CB(contigs, 1) && block_overlap(CB(contigs, 1), &twoblk))
                )){
                missing = true;
            }
        }
        for(int i = 0; i < contigs->size; i++){
            qblk = contigs->block[i];
            if(qblk){
                tcon = QT_SGC(syn, qblk);
                tblk = QT_SGCB(syn, qblk);
                printf("%s %s %u %u %d\n",
                       seqname, tcon->name, tblk->start, tblk->stop, missing);
            }
            else {
                printf("%s %s %s %s %d\n", 
                       seqname, tcon->name, "NA", "NA", missing);
            }
        }
    }
}

void analysis_filter(Synmap * syn, FILE * hitfile,
                     bool(*classifier)(Synmap *, Link *, void *), void * arg){
    Link link;
    char * line;
    size_t len;
    int read;
    bool agrees;
    while ((read = getline(&line, &len, hitfile)) != EOF){
        scanf(line, "%d %d %d %d %d %d\n",
                     &link.qseqid, &link.qstart, &link.qstop,
                     &link.tseqid, &link.tstart, &link.tstop);
        agrees = classifier(syn, &link, arg);
        printf("%d %s\n", agrees, line);
    }
}

bool single_advocate(Synmap * syn, Link * query, void * width_ptr){
    uint width, max_pos, min_pos;
    Contig * con;
    Block * qblk, *tblk;
    width = * (uint *)width_ptr; 

    con = get_region(SGC(syn, 0, query->qseqid),
                                  query->qstart, query->qstop);

    if(!con->start_sorted)
        sort_blocks_by_start(con);

    if(!con->stop_sorted)
        sort_blocks_by_stop(con);

    min_pos = (query->qstart > width) ? query->qstart - width : 0;
    max_pos = query->qstop + width;

    // look down
    int lo_id = CB_STOPID(con, 0);
    for(; lo_id >= 0 && CB_STOP(con, lo_id) > min_pos; lo_id--){
        qblk = CB(con, lo_id); 
        if(qblk->stop >= query->qstart)
            continue;
        if(qblk->oseqid == query->tseqid){
            tblk = QT_SGCB(syn, qblk);
            if(overlap(query->tstart, query->tstop, tblk->start, tblk->stop))
                return true;
        }
    }

    // look up 
    int hi_id = CB_STARTID(con, con->size - 1);
    for(; hi_id < con->size && CB_START(con, hi_id) < max_pos; hi_id++){
        qblk = CB(con, hi_id); 
        if(qblk->start <= query->qstop)
            continue;
        if(qblk->oseqid == query->tseqid){
            tblk = QT_SGCB(syn, qblk);
            if(overlap(query->tstart, query->tstop, tblk->start, tblk->stop))
                return true;
        }
    }
    return false;
}
