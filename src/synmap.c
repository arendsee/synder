#include <assert.h>
#include <errno.h>

#include "synmap.h"

Synmap *init_Synmap()
{
  Synmap *syn = (Synmap *) malloc(sizeof(Synmap));
  syn->genome = (Genome **) malloc(2 * sizeof(Genome *));
  return (syn);
}

void free_Synmap(Synmap * synmap)
{
  if (synmap != NULL) {
    free_Genome(SG(synmap, 0));
    free_Genome(SG(synmap, 1));
    free(synmap->genome);
    free(synmap);
  }
}

void print_Synmap(Synmap * synmap, bool forward)
{
  print_Genome(SG(synmap, 0), forward);
  print_Genome(SG(synmap, 1), forward);
}

void sort_all_contigs(Synmap * synmap)
{
  for (int genid = 0; genid < 2; genid++) {
    for (int conid = 0; conid < SG(synmap, genid)->size; conid++) {
      sort_blocks_by_start(SGC(synmap, genid, conid));
      sort_blocks_by_stop(SGC(synmap, genid, conid));
    }
  }
}

void set_overlap_group(Synmap * syn){
  sort_all_contigs(syn);

  // Holds current overlapping group id
  size_t grpid = 0;
  // Needed for determining overlaps and thus setids
  size_t maximum_stop = 0;
  // The stop position of the current interval
  size_t this_stop = 0;
  // Block temporary holder
  Block * blk;

  // Loop through target and query genomes
  // g := genome id (0 is query, 1 is target)
  for (int g = 0; g <= 1; g++){
    // Loop through each contig in the query genome
    // i := contig id
    for (size_t i = 0; i < SG(syn,g)->size; i++) {
      maximum_stop = 0;
      // Loop though each block in current contig
      // j := block id
      for (size_t j = 0; j < SGC(syn,g,i)->size; j++){
        blk = SGCB(syn,g,i,j);
        this_stop = blk->pos[1];
        // If the start is greater than the maximum stop, then the block is in
        // a new adjacency group. For this to work, Contig->block must be
        // sorted by start. This sort is performed in build_tree.
        if(blk->pos[0] > maximum_stop){
            grpid++;
        }
        if(this_stop > maximum_stop){
          maximum_stop = this_stop;
        }
        blk->grpid = grpid;
      }
      // increment to break adjacency between contigs and genomes
      grpid++;
    }
  }
}


// Link blocks to nearest non-overlapping up and downstream blocks
// For example, given these for blocks:
//  |---a---|
//            |--b--|
//             |----c----|
//                     |---d---|
//                               |---e---|
// a->adj := (NULL, b)
// b->adj := (a, e)
// c->adj := (a, e)
// d->adj := (a, e)
// e->adj := (d, NULL)
void link_adjacent_blocks_directed(Contig * con, Direction d){
  size_t hi_idx = REL_LO_IDX(con, d); // downstream blocks
  size_t lo_idx = hi_idx;             // upstream blocks
  // 0 or (con->size - 1)
  size_t N = REL_HI_IDX(con, d);
  Block ** lo_blks = d ? con->by_stop : con->block;
  Block ** hi_blks = d ? con->block : con->by_stop;
  Block * lo, * hi;
  Block * lo_next = NULL;
  while(true){
    lo = lo_blks[lo_idx];
    hi = hi_blks[hi_idx];
    if(lo_idx != N){
      lo_next = REL_NEXT(lo_blks, lo_idx, d);
    }

    // In diagrams:
    // <--- indicates a hi block
    // ---> indicates a lo block
    // All diagrams are relative to the d==HI direction

    //       ---> 
    // <---
    // OR
    // --->
    //   <---
    // This should should occur only at the beginning
    if(REL_LT(hi->pos[!d], lo->pos[d], d)){
        hi->adj[!d] = NULL;
        if(hi_idx == N)
            break;
        else
            REL_INC(hi_idx, d);
    }
    //  lo     next
    // ---->  ---->
    //               <--- 
    // If next is closer, and overlapping the hi, increment lo
    // So you increment lo until it is adjacent to the current hi
    else if(REL_LT(lo_next->pos[d], hi->pos[!d], d) &&
            lo_next->grpid != hi->grpid)
    {
        if(lo_idx != N)
            REL_INC(lo_idx, d);
    }
    // --->
    //      <---
    // The current lo is next to, and not overlapping, current hi
    else {
        hi->adj[!d] = lo;
        if(hi_idx == N)
            break;
        else
            REL_INC(hi_idx, d);
    }
  }
}
void link_adjacent_blocks(Synmap * syn){
  for (int genid = 0; genid <= 1; genid++){
    for (size_t conid = 0; conid < SG(syn, genid)->size; conid++) {
      link_adjacent_blocks_directed(SGC(syn, genid, conid), HI);
      link_adjacent_blocks_directed(SGC(syn, genid, conid), LO);
    }
  }
}

// ---- A local utility structure used to build contiguous sets ----
typedef struct Node{
  struct Node * down;
  Block * blk; 
} Node;
Node * init_node(Block * blk)
{
  Node * node = (Node *)malloc(sizeof(Node));
  node->blk = blk;
  node->down = NULL;
  return(node);
}
void remove_node(Node * node)
{
  if(node->down != NULL){
    Node * tmp = node->down;
    node->blk = node->down->blk;
    node->down = node->down->down;
    free(tmp);
  }
}
void free_node(Node * node)
{
  if(node->down != NULL)
    free_node(node->down);
  free(node);
}

void link_contiguous_blocks(Synmap * syn)
{
  Block * blk;
  Node * node;
  int qdiff, tdiff;
  char this_strand, older_strand;
  size_t setid = 0;
  for (size_t i = 0; i < SG(syn,0)->size; i++) {
    node = init_node(SGCB(syn, 0, i, 0));
    for (size_t j = 1; j < SGC(syn,0,i)->size; j++){
      blk = SGCB(syn, 0, i, j);
      this_strand = blk->over->strand;
      while(true){
        qdiff = blk->grpid - node->blk->grpid;
        tdiff = blk->over->grpid - node->blk->over->grpid;
        older_strand = node->blk->over->strand;
        // If adjacent
        if((qdiff == 1) &&
           (this_strand == older_strand) &&
           ((tdiff ==  1 && this_strand == '+') ||
            (tdiff == -1 && this_strand == '-')))
        {
          blk->setid = node->blk->setid;
          blk->cnr[0] = node->blk;
          node->blk->cnr[1] = blk;
          node->blk = blk;
          // TODO: in strange cases, one node might be adjacent to multiple
          // nodes. By placing a break here, I just take the first. I need
          // explicit handling for this case.
          break;
        }
        // If at bottom
        else if(node->down == NULL){
          blk->setid = ++setid;
          node->down = init_node(blk);
          break;
        }
        // If definitely not adjacent
        else if(qdiff > 1){
          remove_node(node);
        }
        else {
          node = node->down;
        }
      }
    }
    free_node(node);
  }
}
