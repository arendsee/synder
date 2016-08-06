#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "contiguous.h"
#include "contiguous-map.h"
#include "contig.h"
#include "synmap.h"

// A list of lists to store contiguous sets
typedef struct CSList{
  struct CSList * down;
  struct CSList * over;
  ContiguousNode * cnode; 
} CSList;
void free_CSList(CSList * cslist);
void print_CSList(CSList * cslist, int level);
void add_cnode_CSList(CSList * cslist, ContiguousNode * cnode);
CSList * init_down_CSList(ContiguousNode * cnode);
CSList * init_over_CSList(ContiguousNode * cnode);
CSList * init_CSList();



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

    //Contig * qcon = SGC(syn, 0, chrid);

    CSList * cslist = init_CSList();

    /* Add all contig to a list of lists
     * Each level contains
     *  1. `over` - pointer to a list of nodes, each holding 1 ContiguousNode object.
     *  2. `down` - pointer to the next contiguous set
     */
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
  else if(cslist->over == NULL){
    cslist->over = init_over_CSList(cnode);
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
