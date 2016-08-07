#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "contiguous.h"
#include "contiguous-map.h"
#include "contig.h"
#include "synmap.h"

// A local utility structure used filter and store contiguous sets
typedef struct CSList{
  struct CSList * next;
  ContiguousNode * lo; 
  ContiguousNode * hi; 
  int setid;
} CSList;
CSList * init_empty_CSList();
CSList * init_CSList(ContiguousNode * cnode);
void add_cnode_CSList(CSList * cslist, ContiguousNode * cnode);
void free_CSList(CSList * cslist);

void contiguous_query(Synmap * syn, FILE * intfile, bool pblock)
{

  // count total number of unique block-block pairs for hashmap
  ContiguousMap *cmap = populate_contiguous_map(syn);

  char seqname[128];
  int chrid, start, stop;
  ResultContig * rc;
  Contig *contig;
  size_t length = 1024;
  char *line = (char *) malloc(length * sizeof(char));
  int set_min, set_max;
  while (fgets(line, length, intfile) && !feof(intfile)) {
    if (!sscanf
        (line, "%d %*s %*s %d %d %*s %*c %*s %s\n", &chrid, &start, &stop,
         seqname)) {
      printf("invalid input\n");
      continue;
    }

    Contig *qcon = SGC(syn, 0, chrid);

    rc = get_region(qcon, start, stop);
    contig = rc->contig;

    CSList * cslist = init_empty_CSList();
    CSList * root = cslist;

    // get list of highest and lowest members of each contiguous set
    for(int i = 0; i < contig->size; i++){
        add_cnode_CSList(cslist, cmap->map[contig->block[i]->linkid]); 
    }

    // TODO, need overlaps of inbetween interval flanks

    Contig * tcon;
    int mask, flag;
    int interval[2] = { 0, 0 };
    bool inverted;
    Block *min_blk, *max_blk, *blk;
    for(; cslist != NULL; cslist = cslist->next){

        set_min = get_min(cslist->lo);
        set_max = get_max(cslist->hi);

        tcon = SGC(syn, 1, cslist->lo->feature->oseqid);
        min_blk = cslist->lo->match;
        max_blk = cslist->hi->match;
        inverted = max_blk->strand == '-';
        mask = 0; // reset the mask

        #define AL  cslist->lo->match->start
        #define ALV cslist->lo->match->stop
        #define AR  cslist->hi->match->stop
        #define ARV cslist->hi->match->start

        #define BL  stop > set_max ? \
                    cslist->lo->match->stop : \
                    cslist->lo->prev->match->stop
        #define BLV stop > set_max ? \
                    cslist->lo->match->start : \
                    cslist->lo->prev->match->start
        #define BR  stop < set_min ? \
                    cslist->lo->match->start : \
                    cslist->hi->next->match->start
        #define BRV stop < set_min ? \
                    cslist->hi->match->stop : \
                    cslist->hi->next->match->stop
        #define CLOSEST_L closest_block_below(tcon, cslist->lo->match->start)
        #define UL  blk->stop
        #define ULV blk->start
        #define CLOSEST_R closest_block_above(tcon, cslist->lo->match->stop)
        #define UR  blk->start
        #define URV blk->stop
        #define XL  0
        #define XLV 1000000000 // I need a better solution for this ...
        #define XR  1000000000
        #define XRV 0

        if(start <= min_blk->stop && start >= min_blk->start){
            mask |= ANCHORED_L;
            if(inverted){
                interval[1] = ALV;
            } else {
                interval[0] = AL;
            }
        }
        else if(start < set_min){
            blk = CLOSEST_L;
            if(blk == NULL){
                mask |= EXTREME_L;
                if(inverted){
                    interval[0] = XLV;
                } else {
                    interval[1] = XL;
                }
            } else {
                mask |= UNBOUND_L;
                if(inverted){
                    interval[0] = ULV;
                } else {
                    interval[1] = UL;
                }
            }
        }
        else{
            mask |= BOUND_L;
            if(inverted){
                interval[1] = BLV - 1;
            } else {
                interval[0] = BL + 1;
            }
        }

        if(stop > max_blk->start && stop < max_blk->stop){
            mask |= ANCHORED_R;
            if(inverted){
                interval[0] = ARV;
            } else {
                interval[1] = AR;
            }
        }
        else if(stop > set_max){
            blk = CLOSEST_R;
            if(blk == NULL){
                mask |= EXTREME_R;
                if(inverted){
                    interval[0] = XRV;
                } else {
                    interval[1] = XR;
                }
            } else {
                mask |= UNBOUND_R;
                if(inverted){
                    interval[0] = URV;
                } else {
                    interval[1] = UR;
                }
            }
        }
        else{
            mask |= BOUND_R;
            if(inverted){
                interval[0] = BRV + 1;
            } else {
                interval[1] = BR - 1;
            }
        }
        #undef AL
        #undef ALV
        #undef AR
        #undef ARV
        #undef BL
        #undef BLV
        #undef BR
        #undef BRV
        #undef CLOSEST_L
        #undef UL
        #undef ULV
        #undef CLOSEST_R
        #undef UR
        #undef URV
        #undef XL
        #undef XLV
        #undef XR
        #undef XRV

        switch(mask){
            case F_AA:
                flag = 0;
                break;
            case F_AB:
                flag = 0;
                break;
            case F_BA:
                flag = 0;
                break;
            case F_BB:
                flag = 0;
                break;
            case F_UA:
                flag = 1;
                break;
            case F_XA:
                flag = 1;
                break;
            case F_AU:
                flag = 2;
                break;
            case F_UU:
                flag = 3;
                break;
            case F_UB:
                flag = 4;
                break;
            case F_BU:
                flag = 5;
                break;
            case F_AX:
                flag = 2;
                break;
            case F_UX:
                flag = 7;
                break;
            case F_BX:
                flag = 7;
                break;
            case F_XB:
                flag = 6;
                break;
            case F_XU:
                flag = 6;
                break;
            case F_XX:
                flag = 8;
                break;
            default:
                flag = -1;
                break;
        }

        printf("%s\t%s\t%i\t%i\t%s\t%i\t%i\t.\t%i\n",
            seqname, qcon->name, start, stop, tcon->name, interval[0], interval[1], flag);
    }

    free(rc);
    free_partial_Contig(contig);
    free_CSList(root);
  }
  free(line);
  free_ContiguousMap(cmap);
}


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
  else if(cslist->lo == NULL){ // first entry of empty CSList
    cslist->lo = cnode;
    cslist->hi = cnode;
    cslist->setid = cnode->setid;
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
