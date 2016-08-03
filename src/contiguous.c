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
                    printf(">\t%u\t", interval); \
                  printf("%s\t%s\t%u\t%u\t%s\t%u\t%u\t%c\t%d\n", \
                  seqname, qcon->name, start, stop, \
                  tcon->name, tblk->start, tblk->stop, '.', flag);

#define PRINT_Q printf("Q\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\n", \
                seqname, qcon->name, q_blk->start, q_blk->stop, \
                t_con->name, t_blk->start, t_blk->stop, interval);

#define PRINT_T printf("T\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\n", \
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

  uint32_t interval = 0;
  uint32_t missloc = 0;
  char seqname[128];
  int chrid, start, stop, flag;
  Contig *contigs;
  Contig *tcon;
  Block *qblk;
  uint32_t blkid;
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
    uint32_t region[2] = { 0, 0 };
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
              && offset > 0 ? (uint32_t) offset : 0;
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
              && offset > 0 ? (uint32_t) offset : 0;
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
            tblk->stop = offset >= 0 ? (uint32_t) offset : 0;
          } else {
            flag = SI_LUNBND;
            tblk->stop = t_blk->start;
            offset =
              stop >= q_blk->stop ? stop - q_blk->stop : q_blk->stop - stop;
            offset = t_blk->start - offset;
            tblk->start = offset > 0 ? (uint32_t) offset : 0;
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
            tblk->start = offset > 0 ? (uint32_t) offset : 0;
          } else {
            flag = SI_RUNBND;
            tblk->start = t_blk->stop;
            offset =
              start <=
              q_blk->start ? q_blk->start - start : start - q_blk->start;
            offset += t_blk->stop;
            tblk->stop = offset >= 0 ? (uint32_t) offset : 0;
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
        printf("I\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\n",
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
              printf("I\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\n",
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

ContiguousMap *populate_contiguous_map(Synmap * syn)
{
  // Sum unique blocks in current synmap
  size_t size = 0;
  for (uint32_t i = 0; i < SG(syn,0)->size; i++) {
    size += SGC(syn,0,i)->size;
  }

  // Initialize new ContiguousMap to that size to serve as hashmap
  // when looking up overlapping blocks
  ContiguousMap *cmap = init_ContiguousMap(size);

  // Loop through each contig in the query genome
  // i := contig id
  for (uint32_t i = 0; i < SG(syn,0)->size; i++) {

    uint32_t overlap_bound[2] = { 0, 0 };

    // Loop though each block in current contig
    // j := block id
    for (uint32_t j = 0; j < SGC(syn,0,i)->size; j++) {

      // Set current block up to be a new node
      ContiguousNode *cnode = init_ContiguousNode(syn, i, j);
      // Initialize flags to normal
      cnode->flag = CNF_NORMAL;

      int clinkid = cnode->feature->linkid;
      int coseqid = cnode->feature->oseqid;
      int coblkid = cnode->feature->oblkid;

      // Start of contig is always a default node
      cmap->map[clinkid] = cnode;
      if (j == 0)
        continue;

      // Otherwise we get to figure out if the current block can be added to or
      // breaks the current ContiguousSet
      
      // Get a ContiguousNode object from a Block id
      #define BLKID2NODE(idx) cmap->map[syn->genome[0]->contig[i]->block[idx]->linkid]

      // Test for overlap on query side
      bool q_overlap = false;
      bool t_overlap = false;
      ContiguousNodeFlag prior_flag = BLKID2NODE(j - 1)->flag;

      // k := overlapping block id
      for (int k = overlap_bound[0]; k <= overlap_bound[1]; k++) {
        q_overlap = block_overlap(cnode->feature, SGCB(syn,0,i,k)) || q_overlap;
        t_overlap =
          ((block_overlap(cnode->match, BLKID2NODE(k)->match)) &&
           (coseqid == SGCB(syn,0,i,k)->oseqid)) || t_overlap;
        if (q_overlap || t_overlap)
          continue;
      }
      if (q_overlap) {
        cnode->flag = CNF_QOVER;
      }
      if (t_overlap) {
        cnode->flag = CNF_TOVER;
      }
      if (q_overlap || t_overlap) {
        if (q_overlap && t_overlap)
          cnode->flag = CNF_BOVER;
        overlap_bound[1] = j;
        continue;
      } else {
        overlap_bound[0] = j;
        overlap_bound[1] = j;
      }

      // On different contig from previous 
      if (coseqid != BLKID2NODE(j - 1)->feature->oseqid) {
        continue;
      }

      // Regular contiguous interval
      else if (coblkid == SGCB(syn,0,i,j-1)->oblkid + 1) {
        cnode->flag = CNF_NORMAL;
        // link this node to the previous
        BLKID2NODE(j - 1)->next = cnode;
        cmap->map[clinkid]->prev = BLKID2NODE(j - 1);
      }

      //Twist to left of previous block
      else if (coblkid < SGCB(syn,0,i,j-1)->oblkid) {
        cnode->flag = (prior_flag == CNF_LEFT || prior_flag == CNF_NORMAL) ? CNF_LEFT : CNF_NORMAL;
        BLKID2NODE(j - 1)->next = NULL;
        if (coblkid == SGCB(syn,0,i,j-1)->oblkid - 1) {
          cnode->flag = CNF_INV;
          BLKID2NODE(j - 1)->next = cnode;
          BLKID2NODE(j - 1)->flag = CNF_INV;
          cmap->map[clinkid]->prev = BLKID2NODE(j - 1);
        }
      }

      // Twist to right, possible transposition
      else if (coblkid > SGCB(syn,0,i,j-1)->oblkid) {
        if (j + 1 < SGC(syn,0,i)->size) {
          cnode->flag = SGCB(syn,0,i,j)->oblkid == SGCB(syn,0,i,j+1)->oblkid + 1 ? CNF_INV : CNF_RIGHT;
        } else {
          cnode->flag = CNF_RIGHT;
        }
      }

      // Default case that should never be reached.
      else {                  
        fprintf(stderr, "Something strange happened\n");
        continue;
      }
    }

    #undef BLKID2NODE

  }
  return cmap;

}
