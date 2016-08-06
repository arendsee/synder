#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "contiguous.h"
#include "contig.h"
#include "synmap.h"
#include "contiguous-utils.c"

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

  print_ContiguousMap(cmap);

  free_ContiguousMap(cmap);

//  char seqname[128];
//  int chrid, start, stop;
//  ResultContig * rc;
//  Contig *contigs;
//  size_t length = 1024;
//  char *line = (char *) malloc(length * sizeof(char));
//  while (fgets(line, length, intfile) && !feof(intfile)) {
//    if (!sscanf
//        (line, "%d %*s %*s %d %d %*s %*c %*s %s\n", &chrid, &start, &stop,
//         seqname)) {
//      printf("invalid input\n");
//      continue;
//    }
//    rc = get_region(SGC(syn, 0, chrid), start, stop);
//    contigs = rc->contig;
//
//    //Contig * qcon = SGC(syn, 0, chrid);
//
//    Block * tblk;
//    Block * min_blk;
//    Block * max_blk;
//
//    CSList * cslist = init_CSList();
//
//    /* Add all contigs to a list of lists
//     * Each level contains
//     *  1. `over` - pointer to a list of nodes holding 1 ContiguousNode object.
//     *  2. `down` - pointer to the next contiguous set
//     */
//    for(int i = 0; i < contigs->size; i++){
//        add_cnode_CSList(cslist, cmap->map[contigs->block[i]->linkid]); 
//    }
//
//    CSList * over;
//    CSList * root = cslist;
//    bool inverse;
//    int mask = 0;
//    int flag = 0;
//    ContiguousNode * lol, * lor, * hil, * hir;
//    for(; cslist != NULL; cslist = cslist->down){
//        lol = lor = hil = hir = NULL;
//        mask = 0;
//        for(; over != NULL; over = over->over){
//            
//        }
//
//        if(start < min_blk->stop && start > min_blk->start){
//            mask |= ANCHORED_L;
//        }
//        else if(start > min_blk->stop && min_cnode->prev != NULL){
//            mask |= BOUND_L;
//        }
//        else{
//            mask |= UNBOUND_L;
//        }
//
//        if(stop > max_blk->start && stop < max_blk->stop){
//            mask |= ANCHORED_R;
//        }
//        else if(stop < max_blk->start && max_cnode->next != NULL){
//            mask |= BOUND_R;
//        }
//        else{
//            mask |= UNBOUND_R;
//        }
//
//        #define SNAP_RS
//        #define SNAP_RSV
//        #define SNAP_RE
//        #define SNAP_REV
//        #define SNAP_LS
//        #define SNAP_LSV
//        #define SNAP_LE
//        #define SNAP_LEV
//    
//        switch(mask){
//            case F_AA:
//                flag = F_AA;
//                printf("AA\n");
//                break;
//            case F_AB:
//                flag = F_AB;
//                printf("AB\n");
//                break;
//            case F_BA:
//                flag = F_BA;
//                printf("BA\n");
//                break;
//            case F_BB:
//                flag = F_BB;
//                printf("BB\n");
//                break;
//            case F_AU:
//                flag = F_AU;
//                printf("AU\n");
//                break;
//            case F_UA:
//                flag = F_UA;
//                printf("UA\n");
//                break;
//            case F_UU:
//                flag = F_UU;
//                printf("UU\n");
//                break;
//            case F_BU:
//                flag = F_BU;
//                printf("BU\n");
//                break;
//            case F_UB:
//                flag = F_UB;
//                printf("UB\n");
//                break;
//            default:
//                fprintf(stderr, "Something strange happened in contiguous.c\n");
//                break;
//
//        }
//    }
//    free(rc);
//    free_CSList(root);
//    free_ContiguousMap(cmap);
//  }

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
  cnode->match = QT_SGCB(syn, cnode->feature);
  cnode->qblkid = blkid;
  return cnode;  
}

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
  for(int i = 0; i < SG(syn, 0)->size; i++){
    size_t nblks = SGC(syn, 0, i)->size;
    for(int j = 0; j < nblks; j++){
      ContiguousNode *cnode = init_ContiguousNode(syn, i, j);
      cmap->map[SGCB(syn, 0, i, j)->linkid] = cnode;
    }
  }

  // adjgrp holds ids for each set of overlapping intervals on each genome.  It
  // indexed by linkid. This allows lookup of block adjacency. Two blocks are
  // adjacent if their setids differ by 1 (if on plus strand) or -1 (if on
  // negative strand).
  int ** adjgrp = (int **)malloc(2 * sizeof(int*));
  // Holds current set id
  size_t n[2] = { 0, 0 };
  // Needed for determining overlaps and thus setids
  size_t maximum_stop = 0;
  // The stop position of the current interval
  size_t this_stop = 0;

  // Loop through target and query genomes
  // g := genome id (0 is query, 1 is target)
  Block * blk;
  for (int g = 0; g <= 1; g++){
    adjgrp[g] = (int *)malloc(size * sizeof(int));
    n[g] = 0;
    // Loop through each contig in the query genome
    // i := contig id
    for (size_t i = 0; i < SG(syn,g)->size; i++) {
      // Loop though each block in current contig
      // j := block id
      for (size_t j = 0; j < SGC(syn,g,i)->size; j++){
        blk = SGCB(syn,g,i,j);
        this_stop = blk->stop;
        if(j == 0){
          n[g]++; 
          maximum_stop = 0;
        }
        // If the start is greater than the maximum stop, then the block is in
        // a new adjacency group. For this to work, Contig->block must be
        // sorted by start. This sort is performed in build_tree.
        else if(blk->start > maximum_stop){
            n[g]++;
        }
        if(this_stop > maximum_stop){
          maximum_stop = this_stop;
        }
        adjgrp[g][blk->linkid] = n[g];
      }
    }
  }

  size_t setid = 0;
  cmap->map[0]->setid = setid;
  NodeList * root = init_NodeList(cmap->map[0]);
  NodeList * node;
  ContiguousNode * this_cnode;
  int qdiff, tdiff;
  char this_strand, older_strand;
  for(int i = 1; i < size; i++){
    // reset node to root
    node = root;
    this_cnode = cmap->map[i];
    this_strand = this_cnode->match->strand;
    while(true){
      qdiff = adjgrp[0][i] - adjgrp[0][node->cnode->feature->linkid];
      tdiff = adjgrp[1][i] - adjgrp[1][node->cnode->feature->linkid];
      older_strand = node->cnode->match->strand;
      // If adjacent
      if((qdiff == 1) &&
         (this_strand == older_strand) &&
         ((tdiff ==  1 && this_strand == '+') ||
          (tdiff == -1 && this_strand == '-')))
      {
        this_cnode->setid = setid;
        this_cnode->prev = node->cnode;
        node->cnode->next = this_cnode;
        node->cnode = this_cnode;
        // TODO: in strange cases, one node might be adjacent to multiple
        // nodes. By placing a break here, I just take the first. I need
        // explicit handling for this case.
        break;
      }
      // If definitely not adjacent
      else if(qdiff > 1){
        node = remove_node_NodeList(node);
        continue;
      }
      // IF at bottom
      else if(node->down == NULL || node->cnode == NULL){
        setid++;
        this_cnode->setid = setid;
        init_down_NodeList(node, this_cnode);
        break;
      }
    }
  }

  free_NodeList(root);
  free(adjgrp[0]);
  free(adjgrp[1]);
  free(adjgrp);

  return cmap;
}
