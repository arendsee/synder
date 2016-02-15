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
    bool missing;
    while ((fscanf(intfile,
                   "%d %*s %*s %d %d %*c %*c %*s %s\n",
                   &chrid, &start, &stop, seqname)) != EOF)
    {
        missing = false;
        contigs = get_region(SGC(syn, 0, chrid), start, stop);
        if(contigs->size < 3){
            qblk = contigs->block[0];
            if(!overlap(qblk->start, qblk->stop, start, stop) || contigs->size == 1){
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
