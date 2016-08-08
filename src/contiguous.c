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
  ContiguousNode * bound[2]; 
  int setid;
} CSList;
CSList * init_empty_CSList();
CSList * init_CSList(ContiguousNode * cnode);
void add_cnode_CSList(CSList * cslist, ContiguousNode * cnode);
void free_CSList(CSList * cslist);


typedef struct BoundResult{
  uint bound;
  int mask; 
} BoundResult;

BoundResult * get_si_bound(
  uint q,
  uint set_bounds[2],
  ContiguousNode * blk_bounds[2],
  Direction d,
  bool inverted
);

void contiguous_query(Synmap * syn, FILE * intfile, bool pblock)
{

  // count total number of unique block-block pairs for hashmap
  ContiguousMap *cmap = populate_contiguous_map(syn);

  uint bounds[2];
  uint set_bounds[2]; // Max and min nodes in contiguous set 
  ContiguousNode * blk_bounds[2]; // Bounds for blocks returned from itree

  BoundResult * bound_results[2];

  CSList *cslist, *root;

  int mask, flag;

  size_t length = 1024;
  char *line = (char *) malloc(length * sizeof(char));
  char seqname[128];
  int chrid;
  ResultContig * rc;
  Contig *contig;
  while (fgets(line, length, intfile) && !feof(intfile)) {
  if (!sscanf
    (line, "%d %*s %*s %d %d %*s %*c %*s %s\n", &chrid, &bounds[LO], &bounds[HI],
     seqname)) {
    printf("invalid input\n");
    continue;
  }

  rc = get_region(SGC(syn, 0, chrid), bounds[LO], bounds[HI]);
  contig = rc->contig;

  cslist = init_empty_CSList();
  root = cslist;

  bool inverted;

  // get list of highest and lowest members of each contiguous set
  for(int i = 0; i < contig->size; i++){
    add_cnode_CSList(cslist, cmap->map[contig->block[i]->linkid]); 
  }

  for(; cslist != NULL; cslist = cslist->next){

    set_bounds[LO] = get_set_bound(cslist->bound[LO], LO);
    set_bounds[HI] = get_set_bound(cslist->bound[HI], HI);

    blk_bounds[LO] = cslist->bound[LO];
    blk_bounds[HI] = cslist->bound[HI];

    inverted = blk_bounds[HI]->match->strand == '-';

    bound_results[inverted ^ LO] =
      get_si_bound(bounds[LO], set_bounds, blk_bounds, LO, inverted);
    bound_results[inverted ^ HI] =
      get_si_bound(bounds[HI], set_bounds, blk_bounds, HI, inverted);

    // The HI and LO bounds each have a 4 bit mask describing the overlap
    // state. Since the states are encoded in the odd bits, shifting one
    // mask one bit over and ORing them together, yields a single 8 bit
    // mask fully describing the joint state.
    mask = ( bound_results[LO]->mask << 1 ) | bound_results[HI]->mask;

    // I determine the overlap cases based on this 8 bit mask, Currently I
    // don't use the other bits in the int for anything, but just to be
    // safe, I'll AND them out.
    //
    // A - anchored
    // B - bound
    // U - unbound
    // X - extreme
    // TODO: decide on a stable set of flags to return to the user, for now, setting all to -1
    switch(mask & 0xFF){
      case F_AA:
        flag = -1;
        break;
      case F_AB:
        flag = -1;
        break;
      case F_BA:
        flag = -1;
        break;
      case F_BB:
        flag = -1;
        break;
      case F_UA:
        flag = -1;
        break;
      case F_XA:
        flag = -1;
        break;
      case F_AU:
        flag = -1;
        break;
      case F_UU:
        flag = -1;
        break;
      case F_UB:
        flag = -1;
        break;
      case F_BU:
        flag = -1;
        break;
      case F_AX:
        flag = -1;
        break;
      case F_UX:
        flag = -1;
        break;
      case F_BX:
        flag = -1;
        break;
      case F_XB:
        flag = -1;
        break;
      case F_XU:
        flag = -1;
        break;
      case F_XX:
        flag = -1;
        break;
      default:
        flag = -1;
        break;
    }

    printf("%s\t%s\t%i\t%i\t%s\t%i\t%i\t.\t%i\n",
      seqname,
      blk_bounds[LO]->feature->parent->name,
      bounds[LO], bounds[HI],
      blk_bounds[LO]->match->parent->name,
      bound_results[LO]->bound, bound_results[HI]->bound,
      flag);
  }
  free(rc);
  free_partial_Contig(contig);
  free_CSList(root);
  }
  free(line);
  free_ContiguousMap(cmap);
}

BoundResult * init_BoundResult(uint bound, int mask){
  BoundResult * br = (BoundResult *)malloc(sizeof(BoundResult));
  br->bound = bound;
  br->mask = mask;
  return br;
}

BoundResult * get_si_bound(
  uint q,
  uint set_bounds[2],
  ContiguousNode * blk_bounds[2],
  Direction d,
  bool inverted)
{
  // Invert orientation mapping to target if search interval is inverted
  Direction vd = inverted ? !d : d;

  // +1 or -1 values for snapping above or below a point
  int offset = d ? 1 : -1;
  // Invert these values if the search interval is inverted
  offset = inverted ? -1 * offset : offset;

  int mask = 0;
  uint bound = 444444; // non-zero to ease debugging

  uint pnt_a = blk_bounds[!d]->feature->pos[!d];

  uint pnt_b = blk_bounds[!d]->feature->pos[ d];

  uint pnt_c = blk_bounds[ d]->feature->pos[!d];

  uint pnt_d = blk_bounds[ d]->feature->pos[ d];

  // All diagrams are shown for the d=HI case, take the mirror image fr d=LO.
  //
  // KEY:
  // |x--  --y| - bounds of the contiguous block; start==x, stop==y 
  // a========b - a syntenic block with start == a and stop == b
  //   <---q    - the query interval, with stop == q (start doesn't matter)
  // a==b--c==d - query bounding blocks in the contiguous set
  // [===]      - a non-bounding block in the same contiguous set
  //  ...  F=== - nearest non-adjacent block ***ON THE TARGET***, F=start
  //      ^     - search interval bound
  //
  // Possible snap positions (relative to query)
  //   |x...[===]-----a=======b-----c=======d-----[===]...y|  ...  F===
  //                 ^        ^    ^        ^    ^                ^


  // This may occur when there is only one element in the ContiguousSet
  //          |x-----y|
  //       ...a=======b      
  //         ^
  //   <---q
  // q < x
  if(REL_LT(q, set_bounds[!d], d)){
    bound = blk_bounds[!d]->match->pos[!vd] - offset;
    mask |= UNBOUND;
  }

  //   |x---[===]-------a=======b...y|   ...    F===
  //                   ^
  //              --q
  // q < a
  else if(REL_LT(q, pnt_a, d)){
    bound = blk_bounds[!d]->match->pos[!vd] - offset;
    mask |= BOUND;
  }

  //   |x...a=======b...y|   ...    F===
  //                ^
  //        <---q
  // q < b
  else if(REL_LE(q, pnt_b, d)){
    bound = blk_bounds[!d]->match->pos[vd];
    mask |= ANCHORED;
  }

  //   |x...a=======b-------c=======d...y|  ...  F===
  //                       ^
  //                  <--q
  // q < c && q > b
  //   (q > b test required since blk_bounds[LO] can equal blk_bounds[HI])
  else if(REL_LT(q, pnt_c, d) && REL_GT(q, pnt_b, d))
  {
    bound = blk_bounds[d]->match->pos[!vd] - offset;
    mask |= BOUND;
  }

  //   |x...a=======b-------c=======d...y|  ...  F===
  //                                ^
  //             <--------------q
  // q < d
  else if(REL_LE(q, pnt_d, d)){
    bound = blk_bounds[d]->match->pos[vd];
    mask |= ANCHORED;
  }

  //   |x...a=======b-------c=======d-------[===]...y|  ...  F===
  //                                       ^
  //              <----------------------q
  // q < y, (which implies there is a node after d)
  else if(REL_LE(q, set_bounds[d], d)){
    bound = blk_bounds[d]->adj[d]->match->pos[!vd] - offset;
    mask |= BOUND;
  }

  // If none of the above, the bound beyond anything in the contiguous set
  // In this case, the hi and lo Contiguous nodes will be the same
  else {
    // Get flanking sequence
    Contig * con = blk_bounds[d]->match->parent;
    Block * nearest =
      closest_block(con, blk_bounds[d]->match->pos[vd], vd);
 
    // nearest block on TARGET side exists
    if(nearest != NULL){

      //    |x...--a=======b|
      //    |x...--c=======d|  ...  F===
      //                           ^
      //                  <-----------q
      // q >= F  ***ON QUERY SIDE***
      if(!inverted && REL_GT(q, nearest->over->pos[!d], d)){
        // TODO: This is a weird case, how should I handle it?
        mask |= ANCHORED;
        bound = nearest->pos[!d] - offset;
      }

      //    |x...--a=======b|
      //    |x...--c=======d|  ...  F===
      //                           ^
      //                    <---q
      else {
        mask |= UNBOUND;
        bound = nearest->pos[!vd] - offset;
      }
    }
    //    |x...--a=======b|
    //    |x...--c=======d|  ...  THE_END
    //                    <---q
    // query is further out than ANYTHING in the synteny map
    else {
      mask |= EXTREME;
      // TODO: Is there a better way to handle unbound extremes?
      bound = d ? 999999999 : 0;
    }
  }

  return init_BoundResult(bound, mask);
}


CSList * init_empty_CSList(){
  CSList * cslist = (CSList *)malloc(sizeof(CSList));
  cslist->next = NULL;
  cslist->bound[LO] = NULL;
  cslist->bound[HI] = NULL;
  cslist->setid = -1; // unset
  return(cslist);
}

CSList * init_CSList(ContiguousNode * cnode){
  CSList * cslist = (CSList *)malloc(sizeof(CSList));
  cslist->next = NULL;
  cslist->bound[LO] = cnode;
  cslist->bound[HI] = cnode;
  cslist->setid = cnode->setid;
  return(cslist);
}

void add_cnode_CSList(CSList * cslist, ContiguousNode * cnode){
  if(cslist->setid == cnode->setid){
    if(cslist->bound[HI] == NULL || cnode->feature->pos[LO] > cslist->bound[HI]->feature->pos[HI]){
        cslist->bound[HI] = cnode;
    }
    if(cslist->bound[LO] == NULL || cnode->feature->pos[HI] < cslist->bound[LO]->feature->pos[LO]){
        cslist->bound[LO] = cnode;
    }
  }
  else if(cslist->bound[LO] == NULL){ // first entry of empty CSList
    cslist->bound[LO] = cnode;
    cslist->bound[HI] = cnode;
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
