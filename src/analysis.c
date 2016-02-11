#include <stdio.h>

#include "synmap.h"
#include "contig.h"
#include "block.h"

void analysis_count(Synmap * syn, FILE * intfile){
    char seqname[128];
    uint count;
    int chrid, start, stop;
    while ((fscanf(intfile,
                   "%d %*s %*s %d %d %*c %*c %*s %s\n",
                   &chrid, &start, &stop, seqname)) != EOF)
    {
        count = count_overlaps(syn->genome[0]->contig[chrid], start, stop);
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
    while ((fscanf(intfile,
                   "%d %*s %*s %d %d %*c %*c %*s %s\n",
                   &chrid, &start, &stop, seqname)) != EOF)
    {
        contigs = get_overlapping(syn->genome[0]->contig[chrid], start, stop);
        for(int i = 0; i < contigs->size; i++){
            qblk = contigs->block[i];
            tcon = syn->genome[1]->contig[qblk->oseqid];
            tblk = tcon->block[qblk->oblkid];
            printf("%s %s %u %u\n", seqname, tcon->name, tblk->start, tblk->stop);
        }
    }
}
