#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "contiguous.h"
#include "contig.h"
#include "synmap.h"

// strand is hardcoded as '.' (strand unknown) because correct strand handling
// is not yet implemented.
#define PRINT_SRC if(pblock) \
                    printf(">\t%lu\t", interval); \
                  printf("%s\t%s\t%u\t%u\t%s\t%u\t%u\t%c\t%d\n", \
                  seqname, qcon->name, start, stop, \
                  tcon->name, tblk->start, tblk->stop, '.', flag);

#define PRINT_Q printf("Q\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%lu\n", \
                seqname, qcon->name, q_blk->start, q_blk->stop, \
                t_con->name, t_blk->start, t_blk->stop, interval);

#define PRINT_T printf("T\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%lu\n", \
               seqname, qcon->name, q_blk->start, q_blk->stop, \
               t_con->name, t_blk->start, t_blk->stop, interval);


ContiguousMap *init_ContiguousMap(size_t size)
{
  ContiguousMap *cmap = (ContiguousMap *) malloc(sizeof(ContiguousMap));
  cmap->size = size;
  cmap->map = (ContiguousNode **) malloc(size * sizeof(ContiguousNode *));

  return cmap;
}

void contiguous_query(Synmap * syn, FILE * intfile, bool pblock)
{
  // count total number of unique block-block pairs for hashmap
  ContiguousMap *cmap = populate_contiguous_map(syn);

  size_t interval = 0;
  size_t missloc = 0;
  char seqname[128];
  int chrid, start, stop, flag;
  Contig *contigs;
  Contig *tcon;
  Block *qblk;
  size_t blkid;
  ContiguousNode *qnode;
  ContiguousNode *original;
  Block twoblk;
  bool missing;
  size_t length = 1024;
  char *line = (char *) malloc(length * sizeof(char));

  while (fgets(line, length, intfile) && !feof(intfile)) {
    if (!sscanf
        (line, "%d %*s %*s %d %d %*s %*c %*s %s\n", &chrid, &start, &stop,
         seqname)) {
      printf("invalid input\n");
      continue;
    }
    contigs = get_region(SGC(syn, 0, chrid), start, stop);
    size_t region[2] = { 0, 0 };
    missing = false;

    // Determine if the results from the region search is an array of valid blocks
    // or is between two blocks
    if (contigs->size == 2) {
      twoblk.start = start;
      twoblk.stop = stop;
      if (!((CB(contigs, 0) && block_overlap(CB(contigs, 0), &twoblk)) ||
            (CB(contigs, 1) && block_overlap(CB(contigs, 1), &twoblk))
          )) {
        if (contigs->block[0] != NULL) {
          missloc = cmap->map[contigs->block[0]->linkid]->qblkid;
        } else {
          missloc = cmap->map[contigs->block[1]->linkid]->qblkid;
          // TODO Q - When will this be needed? Can qblkid be negative? Or is
          // this to deal with overflows?
          missloc = missloc > 0 ? missloc : 0;
        }
        missing = true;
      }
    }

    // Set range of blocks query overlaps. Doing it this way adverts the chance
    // that get_region doesn't return every overlapping block.
    for (int i = 0; i < contigs->size; i++) {
      qblk = contigs->block[i];
      if (qblk != NULL) {
        blkid = cmap->map[qblk->linkid]->qblkid;
        if (i == 0) {
          region[0] = blkid;
          region[1] = blkid;
        } else {
          if (blkid < region[0]) {
            region[0] = blkid;
          }
          if (blkid > region[1]) {
            region[1] = blkid;
          }
        }
      }
    }

    // Loop through the query region bounds and print
    // overlapping target regions
    for (int i = region[0]; i <= region[1]; i++) {

      flag = SI_GOOD;
      qblk = SGCB(syn, 0, chrid, i);
      Contig *qcon = SGC(syn, 0, chrid);
      Block *tblk = init_Block(QT_SGCB(syn, qblk)->start, QT_SGCB(syn, qblk)->stop,
                               0, 0, 0, '.');
      tcon = QT_SGC(syn, qblk);
      Block *q_blk;
      Block *t_blk;
      Contig *t_con;

      // Case "E" examples
      if (missing) {
        if (pblock) {
          q_blk = SGCB(syn, 0, chrid, missloc);
          t_blk = QT_SGCB(syn, q_blk);
          t_con = QT_SGC(syn, q_blk);
          PRINT_Q
          q_blk = SGCB(syn, 0, chrid, missloc + 1);
          t_blk = QT_SGCB(syn, q_blk);
          t_con = QT_SGC(syn, q_blk);
          PRINT_Q
        }
        q_blk = SGCB(syn, 0, chrid, missloc);
        t_blk = QT_SGCB(syn, q_blk);
        tcon = QT_SGC(syn, q_blk);
        int64_t offset;
        if (start < q_blk->start) {     // query region is before block
          if (cmap->map[q_blk->linkid]->flag > CNF_LEFT) {
            // return from start of block, to offest to start of query
            // on target side
            flag = SI_LUNBND;
            tblk->stop = t_blk->start;
            offset = t_blk->start - (q_blk->start - start);
            tblk->start = offset < tblk->stop
              && offset > 0 ? (size_t) offset : 0;
          } else {
            flag = SI_RUNBND;
            tblk->start = t_blk->stop;
            tblk->stop = t_blk->stop + (q_blk->start - start);
          }
        } else {                // query region after block
          if (cmap->map[q_blk->linkid]->flag > CNF_LEFT) {
            // return from end of block, to offest to end of query
            // on target side
            flag = SI_RUNBND;
            tblk->start = t_blk->stop;
            tblk->stop = t_blk->stop + (stop - q_blk->stop);
          } else {
            flag = SI_LUNBND;
            tblk->stop = t_blk->start;
            offset = t_blk->start - (stop - q_blk->stop);
            tblk->start = offset < tblk->stop
              && offset > 0 ? (size_t) offset : 0;
          }
        }
        PRINT_SRC

        blkid = cmap->map[q_blk->linkid]->qblkid;
        if (start < q_blk->start && blkid > 0) {
          q_blk = SGCB(syn, 0, chrid, blkid - 1);
          t_blk = QT_SGCB(syn, q_blk);
          tcon = QT_SGC(syn, q_blk);
          if (cmap->map[q_blk->linkid]->flag > CNF_LEFT) {
            flag = SI_RUNBND;
            tblk->start = t_blk->stop;
            offset =
              stop >= q_blk->stop ? stop - q_blk->stop : q_blk->stop - stop;
            offset += t_blk->stop;
            tblk->stop = offset >= 0 ? (size_t) offset : 0;
          } else {
            flag = SI_LUNBND;
            tblk->stop = t_blk->start;
            offset =
              stop >= q_blk->stop ? stop - q_blk->stop : q_blk->stop - stop;
            offset = t_blk->start - offset;
            tblk->start = offset > 0 ? (size_t) offset : 0;
          }
          PRINT_SRC

        } else if (blkid + 1 < qcon->size) {    // query region after block
          q_blk = SGCB(syn, 0, chrid, blkid + 1);
          t_blk = QT_SGCB(syn, q_blk);
          tcon = QT_SGC(syn, q_blk);
          if (cmap->map[q_blk->linkid]->flag > CNF_LEFT) {
            flag = SI_LUNBND;
            tblk->stop = t_blk->start;
            offset =
              start <=
              q_blk->start ? q_blk->start - start : start - q_blk->start;
            offset = t_blk->start - offset;
            tblk->start = offset > 0 ? (size_t) offset : 0;
          } else {
            flag = SI_RUNBND;
            tblk->start = t_blk->stop;
            offset =
              start <=
              q_blk->start ? q_blk->start - start : start - q_blk->start;
            offset += t_blk->stop;
            tblk->stop = offset >= 0 ? (size_t) offset : 0;
          }
          PRINT_SRC

        }
        // query region is before block

        interval++;
        free_Block(tblk);
        break;
      }

      qnode = cmap->map[qblk->linkid];
      original = cmap->map[qblk->linkid];

      // Set return region assumes case D;
      tblk->start = qnode->match->start;
      tblk->stop = qnode->match->stop;

      // Start is BEFORE current query Block;
      if (start < qblk->start) {
        //Move down contiguous block, stopping at leftmost possible point
        while (qnode->prev != NULL && start > qnode->feature->stop) {
          qnode = qnode->prev;
        }
        if (qnode->flag > CNF_NORMAL && start > qnode->feature->stop) {  //avoid duplicates in cases of overlap
          continue;
        }
        if (start < qnode->feature->start) {    // Check we didn't advance into a case E,F Situation
          if (qnode->flag > CNF_LEFT
              || (original->next != NULL && original->next->flag == CNF_NORMAL)) {
            // The next check here and in the next block is to detect edge cases where block we are checking is
            // registering as a left translation, but is part of a continuous run to the right
            // most often occcurs in the middle of messy overlap blocks
            flag = SI_LUNDET;
            tblk->start = qnode->match->start;
          } else {
            flag = SI_RUNDET;
            tblk->stop = qnode->match->stop;
          }
        } else if (start > qnode->feature->start) {     //Start is contained within current block C,D
          if (qnode->flag > CNF_LEFT) {
            tblk->start = qnode->match->start;
          } else {
            tblk->stop = qnode->match->stop;
          }
        } else {                //Case A,B situations
          if (qnode->flag > CNF_LEFT) {
            tblk->start = qnode->match->stop;
          } else {
            tblk->stop = qnode->match->start;
          }
        }
      }

      if (pblock) {
        //qnode = cmap->map[qblk->linkid];
        q_blk = SGCB(syn, 0, chrid, qnode->qblkid > 0 ? qnode->qblkid - 1 : 0);
        t_blk = QT_SGCB(syn, q_blk);
        t_con = QT_SGC(syn, q_blk);
        PRINT_Q
        t_blk = SGCB(syn, 1, qnode->feature->oseqid,
                     (qnode->feature->oblkid >
                      0 ? qnode->feature->oblkid - 1 : 0));
        q_blk = SGCB(syn, 0, t_blk->oseqid, t_blk->oblkid);
        t_con = SGC(syn, 1, q_blk->oseqid);
        PRINT_T
        printf("I\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%lu\n",
               seqname, qcon->name, qblk->start, qblk->stop,
               tcon->name, qnode->match->start, qnode->match->stop, interval);
      }

      //Stop is AFTER current query Block
      if (stop > qblk->stop) {  //Stop is AFTER
        while (qnode->next != NULL && stop) {
          ContiguousNode *temp = qnode->next;
          //Advance i to avoid repeated ranges due to continuity.
          if (stop > temp->feature->start) {
            i++;
            qnode = temp;
            if (pblock) {
              t_con = SGC(syn, 1, qnode->feature->oseqid);
              printf("I\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%lu\n",
                     seqname, qcon->name, qnode->feature->start,
                     qnode->feature->stop, t_con->name, qnode->match->start,
                     qnode->match->stop, interval);
            }
          } else {
            break;
          }
        }
        if (qnode->flag > CNF_NORMAL && stop < qnode->feature->start) {  //avoid duplicates in cases of overlap
          continue;
        }
        if (stop > qnode->feature->stop) {
          //Case E,F situations
          if (qnode->next == NULL) {
            if (qnode->flag > CNF_LEFT) {
              flag = flag == SI_LUNDET ? SI_BUNDET : SI_RUNDET;
              tblk->stop = qnode->match->stop;
            } else {
              flag = flag == SI_RUNDET ? SI_BUNDET : SI_LUNDET;
              tblk->start = qnode->match->start;
            }
          }
          //Case A,B situations
          else {
            q_blk = SGCB(syn, 0, chrid, i + 1 < qcon->size ? i + 1 : i);
            t_blk = QT_SGCB(syn, q_blk);
            if (cmap->map[qblk->linkid]->flag > CNF_LEFT || qnode->next->flag == CNF_NORMAL) { //similar check to above 
              tblk->stop = t_blk->start;
            } else {
              tblk->start = t_blk->stop;
            }
          }
        }
        //Case C,D 
        else {
          tblk->stop = qnode->match->stop;
          if (qnode->flag > CNF_LEFT) {
            tblk->stop = qnode->match->stop;
          } else {
            tblk->start = qnode->match->stop;
          }
        }
      }

      if (pblock) {
        if (stop < qnode->feature->start) {
          q_blk = qnode->feature;
        } else {
          q_blk = SGCB(syn, 0, chrid, i + 1 < qcon->size ? i + 1 : i);
        }
        t_blk = QT_SGCB(syn, q_blk);
        t_con = QT_SGC(syn, q_blk);
        PRINT_Q
        t_blk =
          SGCB(syn, 1, qnode->feature->oseqid,
               qnode->feature->oblkid + 1 < tcon->size ?
               qnode->feature->oblkid + 1 :
               qnode->feature->oblkid);
        q_blk = SGCB(syn, 0, t_blk->oseqid, t_blk->oblkid);
        t_con = SGC(syn, 1, q_blk->oseqid);
        PRINT_T
      }

      PRINT_SRC

      interval++;

      free_Block(tblk);
    }
    free(contigs->name);
    free(contigs->block);
    free(contigs);
  }
  free(line);
  free_ContiguousMap(cmap);
}

void free_ContiguousMap(ContiguousMap * cmap)
{
  for (int i = 0; i < cmap->size; i++) {
    if (cmap->map[i])
      free(cmap->map[i]);
  }
  if (cmap->map)
    free(cmap->map);
  if (cmap != NULL)
    free(cmap);
}

ContiguousNode * init_ContiguousNode(Synmap * syn, size_t conid, size_t blkid){
  ContiguousNode *cnode = (ContiguousNode *) malloc(sizeof(ContiguousNode));
  cnode->feature = syn->genome[0]->contig[conid]->block[blkid];
  cnode->flag = CNF_UNSET;
  cnode->next = NULL;
  cnode->prev = NULL;
  cnode->match = SGCB(syn, 1, cnode->feature->oseqid, cnode->feature->oblkid);
  cnode->qblkid = blkid;
  return cnode;  
}

// ----------------------------------------------------------------------
// ---- A local utility data structure used to build contiguous sets ----
// Hold list of open ContiguousNode as they are added
typedef struct Wheel {
  ContiguousNode * cnode;
  struct Wheel * next;
  struct Wheel * prev;
} Wheel;
// Reinvent the wheel
Wheel * init_wheel(){
  Wheel * wheel = (Wheel *)malloc(sizeof(Wheel));
  wheel->cnode = NULL;
  wheel->next = wheel;
  wheel->prev = wheel;
  return(wheel);
}
// Remove from Wheel list after the contiguous set is free
Wheel * remove_wheel(Wheel * wheel){
  if(wheel->next == wheel){
    wheel->cnode = NULL;
  } else {
    Wheel * to_free = wheel;
    wheel->prev->next = wheel->next;
    wheel->next->prev = wheel->prev;
    wheel = wheel->next;
    free(to_free);
  }
  return wheel;
}
// Add new ContiguousNode to contiguous set
Wheel * add_cnode(Wheel * wheel, ContiguousNode * cnode){
  wheel->cnode->next = cnode;
  cnode->prev = wheel->cnode;
  wheel->cnode = cnode;
  return wheel->next;
}
// Initialize a new contiguous set
Wheel * add_wheel(Wheel * wheel, ContiguousNode * cnode){
  if(wheel->cnode == NULL){
    wheel->cnode = cnode;
    wheel->next = wheel; // should be unnecessary
    wheel->prev = wheel; // should be unnecessary
  } else {
    Wheel * new_wheel = (Wheel *)malloc(sizeof(Wheel));
    new_wheel->cnode = cnode;
    new_wheel->prev = wheel;
    wheel->next = new_wheel;
    wheel->prev = new_wheel;
  }
  return wheel;
}
// ----------------------------------------------------------------------

ContiguousMap *populate_contiguous_map(Synmap * syn)
{

  // Sum unique blocks in current synmap
  size_t size = 0;
  for (size_t i = 0; i < SG(syn,0)->size; i++) {
    size += SGC(syn,0,i)->size;
  }
  
  // Initialize new ContiguousMap to that size to serve as hashmap
  // when looking up overlapping blocks
  ContiguousMap *cmap = init_ContiguousMap(size);
  
  size_t ** adjgrp = (size_t **)malloc(2 * sizeof(size_t*));
  size_t n[2] = { 0, 0 };
  size_t maximum_stop = 0;
  size_t this_stop = 0;
  
  // Loop through target and query genomes
  // g := genome id (0 is query, 1 is target)
  for (int g = 0; g <= 1; g++){
    adjgrp[g] = (size_t *)malloc(size * sizeof(size_t));
    n[g] = 0;
    // Loop through each contig in the query genome
    // i := contig id
    for (size_t i = 0; i < SG(syn,g)->size; i++) {
      // Loop though each block in current contig
      // j := block id
      for (size_t j = 0; j < SGC(syn,g,i)->size; j++){
        ContiguousNode *cnode = init_ContiguousNode(syn, i, j);
        this_stop = SGCB(syn,g,i,j)->stop;
        if(j == 0){
          n[g]++; 
          maximum_stop = 0;
        }
        // If the start is greater than the maximum stop, then the block is in
        // a new adjacency group. For this to work, Contig->block must be sorted
        // by start. This sort is performed in build_tree.
        else if(SGCB(syn,g,i,j)->start > maximum_stop){
            n[g]++;
        }
        if(this_stop > maximum_stop){
          maximum_stop = this_stop;
        }
        adjgrp[g][SGCB(syn,g,i,j)->linkid] = n[g];
        cmap->map[SGCB(syn,g,i,j)->linkid] = cnode;
      }
    }
  }

  Wheel * wheel = init_wheel();
  Wheel * initial;
  size_t qdiff, tdiff;
  char this_strand;
  bool same_strand;
  bool has_match;
  for(int i = 1; i < size; i++){
    initial = wheel;
    this_strand = QT_SGCB(syn, cmap->map[i]->feature)->strand;
    has_match = false;
    do {
      tdiff = adjgrp[1][i] - adjgrp[1][wheel->cnode->feature->linkid];
      qdiff = adjgrp[0][i] - adjgrp[0][wheel->cnode->feature->linkid];
      same_strand = this_strand == QT_SGCB(syn, wheel->cnode->feature)->strand;
      // If the blocks are adjacent
      if(qdiff == 1 &&
         ((tdiff == 1 && same_strand) || (tdiff == -1 && !same_strand))
      ){
        wheel = add_cnode(wheel, cmap->map[i]);
        has_match = true;
      }
      // If the blocks are not adjacent AND no subsequent ones can be
      else if(qdiff > 1){
        initial = wheel == initial ? wheel->prev : initial;
        wheel = remove_wheel(wheel);
      }
    // Stop if return to initial position
    } while(wheel != initial);
    if(!has_match){
      wheel = add_wheel(wheel, cmap->map[i]);
    }
  }

  return cmap;
}
