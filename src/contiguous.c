#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "contiguous.h"
#include "contig.h"
#include "synmap.h"

#define PRINT_SRC \
  printf("%s\t%s\t%u\t%u\t%s\t%u\t%u\t%c\t%d\n", \
  seqname, qcon->name, start, stop, \
  tcon->name, tblk->start, tblk->stop, '.', flag);


ContiguousMap *init_ContiguousMap(size_t size)
{
  ContiguousMap *cmap = (ContiguousMap *) malloc(sizeof(ContiguousMap));
  cmap->size = size;
  cmap->map = (ContiguousNode **) malloc(size * sizeof(ContiguousNode *));
  return cmap;
}

// ----------------------------------------------------------------------
// ---- A local utility data structure to hold search sets --------------
typedef struct CSList{
  struct CSList * down;
  struct CSList * over;
  ContiguousNode * cnode; 
} CSList;

CSList * init_CSList(){
  CSList * cslist = (CSList *)malloc(sizeof(CSList));
  cslist->cnode = NULL;
  cslist->down = NULL;
  cslist->over = NULL;
  return(cslist);
}

CSList * init_over_CSList(ContiguousNode * cnode){
  CSList * cslist = init_CSList();
  cslist->cnode = cnode;
  return(cslist);
}

CSList * init_down_CSList(ContiguousNode * cnode){
  CSList * cslist = init_CSList();
  cslist->over = init_over_CSList(cnode);
  return(cslist);
}

void add_cnode_CSList(CSList * cslist, ContiguousNode * cnode){
  if(cslist->over != NULL && cslist->over->cnode->setid == cnode->setid){
    CSList * newlist = init_over_CSList(cnode);
    newlist->over = cslist->over;
    cslist->over = newlist;
  }
  else if(cslist->down == NULL){
    cslist->down = init_down_CSList(cnode);
  }
  else{
    add_cnode_CSList(cslist->down, cnode);
  }
}

void free_CSList(CSList * cslist){
  if(cslist->over != NULL){
    free_CSList(cslist->over);
  }
  if(cslist->down != NULL){
    free_CSList(cslist->down);
  }
  free(cslist);
}

void print_CSList(CSList * cslist, int level){
  if(cslist->over != NULL){
    print_CSList(cslist->over, level);
  }
  if(cslist->cnode!= NULL){
    printf(
      "%i (%u, %u)\n",
      level,
      cslist->cnode->feature->start,
      cslist->cnode->feature->stop
    );
  }
  if(cslist->down != NULL){
    print_CSList(cslist->down, level++);
  }
}
// ----------------------------------------------------------------------

void contiguous_query(Synmap * syn, FILE * intfile, bool pblock)
{
  // count total number of unique block-block pairs for hashmap
  ContiguousMap *cmap = populate_contiguous_map(syn);

  char seqname[128];
  int chrid, start, stop;
  ResultContig * rc;
  Contig *contigs;
  size_t length = 1024;
  char *line = (char *) malloc(length * sizeof(char));

  while (fgets(line, length, intfile) && !feof(intfile)) {
    if (!sscanf
        (line, "%d %*s %*s %d %d %*s %*c %*s %s\n", &chrid, &start, &stop,
         seqname)) {
      printf("invalid input\n");
      continue;
    }
    rc = get_region(SGC(syn, 0, chrid), start, stop);
    contigs = rc->contig;

    //Contig * qcon = SGC(syn, 0, chrid);

    CSList * cslist = init_CSList();

    for(int i = 0; i < contigs->size; i++){
        add_cnode_CSList(cslist, cmap->map[contigs->block[i]->linkid]); 
    }

    print_CSList(cslist, 0);

    free(rc);

    free_CSList(cslist);

    free_ContiguousMap(cmap);
  }
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
typedef struct Wheel
{
  ContiguousNode * cnode;
  struct Wheel * next;
  struct Wheel * prev;
} Wheel;
// Reinvent the wheel
Wheel * init_wheel()
{
  Wheel * wheel = (Wheel *)malloc(sizeof(Wheel));
  wheel->cnode = NULL;
  wheel->next = wheel;
  wheel->prev = wheel;
  return(wheel);
}
// Remove from Wheel list after the contiguous set is free
Wheel * remove_wheel(Wheel * wheel)
{
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
Wheel * add_cnode(Wheel * wheel, ContiguousNode * cnode)
{
  wheel->cnode->next = cnode;
  cnode->prev = wheel->cnode;
  wheel->cnode = cnode;
  return wheel->next;
}
// Initialize a new contiguous set
Wheel * add_wheel(Wheel * wheel, ContiguousNode * cnode, size_t setid)
{
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
  wheel->cnode->setid = setid;
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
  size_t setid = 0;
  char this_strand;
  bool same_strand;
  bool has_match;
  for(int i = 0; i < size; i++){
    initial = wheel;
    this_strand = QT_SGCB(syn, cmap->map[i]->feature)->strand;
    has_match = false;
    while(wheel != initial){
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
    }
    if(!has_match){
      wheel = add_wheel(wheel, cmap->map[i], setid);
      setid++;
    }
  }

  return cmap;
}
