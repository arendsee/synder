#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "contiguous.h"
#include "contiguous-map.h"
#include "contig.h"
#include "synmap.h"

// ---- A local utility structure used filter and store contiguous sets ----
typedef struct CSList{
  struct CSList * next;
  ContiguousNode * lo; 
  ContiguousNode * hi; 
  int setid;
} CSList;

CSList * init_empty_CSList(){
  CSList * cslist = (CSList *)malloc(sizeof(CSList));
  cslist->next = NULL;
  cslist->hi = NULL;
  cslist->lo = NULL;
  cslist->setid = -1; // unset
  return(cslist);
}

CSList * init_CSList(ContiguousNode * cnode){
  CSList * cslist = (CSList *)malloc(sizeof(CSList));
  cslist->next = NULL;
  cslist->hi = cnode;
  cslist->lo = cnode;
  cslist->setid = cnode->setid;
  return(cslist);
}

void add_cnode_CSList(CSList * cslist, ContiguousNode * cnode){
  if(cslist->setid == cnode->setid){
    if(cslist->hi == NULL || cnode->feature->start > cslist->hi->feature->start){
        cslist->hi = cnode;
    }
    if(cslist->lo == NULL || cnode->feature->stop < cslist->lo->feature->stop){
        cslist->lo = cnode;
    }
  }
  else if(cslist->next == NULL){
    cslist->next = init_CSList(cnode);
  }
  else{
    add_cnode_CSList(cslist->next, cnode);
  }
}

void free_CSList(CSList * cslist){
  if(cslist->next != NULL)
    free_CSList(cslist->next);
  free(cslist);
}
// -----------------------------------------------------------



void contiguous_query(Synmap * syn, FILE * intfile, bool pblock)
{

  // count total number of unique block-block pairs for hashmap
  ContiguousMap *cmap = populate_contiguous_map(syn);

  print_ContiguousMap(cmap);

  char seqname[128];
  int chrid, start, stop;
  ResultContig * rc;
  Contig *contig;
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
    contig = rc->contig;

    CSList * cslist = init_empty_CSList();

    // get list of highest and lowest members of each contiguous set
    for(int i = 0; i < contig->size; i++){
        add_cnode_CSList(cslist, cmap->map[contig->block[i]->linkid]); 
    }

    free(rc);
    free_partial_Contig(contig);
    free_CSList(cslist);
  }
  free(line);
  free_ContiguousMap(cmap);

}

/*
    Block * tblk;
    Block * min_blk;
    Block * max_blk;
 
    CSList * over;
    bool inverse;
    int mask = 0;
    int flag = 0;
    ContiguousNode * lol, * lor, * hil, * hir;
    for(; cslist != NULL; cslist = cslist->down){
        lol = lor = hil = hir = NULL;
        mask = 0;
        for(; over != NULL; over = over->over){
 
        }
 
        if(start < min_blk->stop && start > min_blk->start){
            mask |= ANCHORED_L;
        }
        else if(start > min_blk->stop && min_cnode->prev != NULL){
            mask |= BOUND_L;
        }
        else{
            mask |= UNBOUND_L;
        }
 
        if(stop > max_blk->start && stop < max_blk->stop){
            mask |= ANCHORED_R;
        }
        else if(stop < max_blk->start && max_cnode->next != NULL){
            mask |= BOUND_R;
        }
        else{
            mask |= UNBOUND_R;
        }
 
        #define SNAP_RS
        #define SNAP_RSV
        #define SNAP_RE
        #define SNAP_REV
        #define SNAP_LS
        #define SNAP_LSV
        #define SNAP_LE
        #define SNAP_LEV
 
        switch(mask){
            case F_AA:
                flag = F_AA;
                printf("AA\n");
                break;
            case F_AB:
                flag = F_AB;
                printf("AB\n");
                break;
            case F_BA:
                flag = F_BA;
                printf("BA\n");
                break;
            case F_BB:
                flag = F_BB;
                printf("BB\n");
                break;
            case F_AU:
                flag = F_AU;
                printf("AU\n");
                break;
            case F_UA:
                flag = F_UA;
                printf("UA\n");
                break;
            case F_UU:
                flag = F_UU;
                printf("UU\n");
                break;
            case F_BU:
                flag = F_BU;
                printf("BU\n");
                break;
            case F_UB:
                flag = F_UB;
                printf("UB\n");
                break;
            default:
                fprintf(stderr, "Something strange happened in contiguous.c\n");
                break;
 
        }
    }
*/
