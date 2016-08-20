#include <assert.h>
#include <errno.h>

#include "synmap.h"

Synmap *init_Synmap()
{
  Synmap *syn = (Synmap *) malloc(sizeof(Synmap));
  syn->size = 2;
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

void link_four_corners__set_head_and_tail(Synmap * syn){
  Contig * con;
  size_t N;
  for (size_t g = 0; g <= 1; g++){
    for (size_t c = 0; c < SG(syn,g)->size; c++) {
      con = SGC(syn, g, c); 
      N = con->size;
      // sort by stop
      sort_blocks(con, true);
      con->head[1] = &con->block[0];
      con->tail[1] = &con->block[N-1];
      for(size_t i = 0; i < N; i++){
        con->block[i].cor[PREV_STOP] = (i == 0)     ? NULL : &con->block[i-1];
        con->block[i].cor[NEXT_STOP] = (i == N - 1) ? NULL : &con->block[i+1];
      }
      // sort by start
      sort_blocks(con, false);
      con->head[0] = &con->block[0];
      con->tail[0] = &con->block[N-1];
      for(size_t i = 0; i < N; i++){
        con->block[i].cor[PREV_START] = (i == 0)     ? NULL : &con->block[i-1];
        con->block[i].cor[NEXT_START] = (i == N - 1) ? NULL : &con->block[i+1];
      }
    }
  }
}

void set_overlap_group(Synmap * syn){

  // Holds current overlapping group id
  size_t grpid = 1;
  // Needed for determining overlaps and thus setids
  long maximum_stop = 0;
  // The stop position of the current interval
  long this_stop = 0;
  // Current Block in linked list
  Block * blk;

  // Loop through target and query genomes
  // g := genome id (0 is query, 1 is target)
  for(size_t g = 0; g <= 1; g++){
    // Loop through each contig in the query genome
    // i := contig id
    for(size_t i = 0; i < SG(syn,g)->size; i++) {
      maximum_stop = 0;
      // Loop through each Block in the linked list
      blk = SGC(syn,g,i)->head[0];
      for(; blk != NULL; blk = blk->cor[1]){
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
  // In diagrams:
  // <--- indicates a hi block
  // ---> indicates a lo block
  // All diagrams and comments relative to the d==HI direction

  Block * lo = con->head[ d]; // first elements by stop
  Block * hi = con->head[!d]; // first elements by start

  // Transformed indices for Block->cor
  // ----------------------------------
  // a   b          c   d
  // <--->   <--->  <--->
  //   a - previous element by start
  //   b - previous element by stop
  //   c - next element by start
  //   d - next element by stop
  //int idx_a = (!d * 2) + !d  // - 0
  //int idx_b = ( d * 2) + !d  // - 2
  int idx_c = (!d * 2) +  d;  // - 1
  int idx_d = ( d * 2) +  d;  // - 3
  
  while(hi != NULL){

    //       ---> 
    // <---
    // OR
    // --->
    //   <---
    // This should should occur only at the beginning
    if(hi->pos[!d] <= lo->pos[d])
    {
      hi->adj[!d] = NULL;
      hi = hi->cor[idx_c];
    }

    //  lo     next
    // ---->  ---->
    //               <---
    // If next is closer, and not overlapping the hi, increment lo
    // You increment lo until it is adjacent to the current hi
    else if(lo->cor[idx_d]->pos[d] < hi->cor[idx_c]->pos[!d])
    {
      lo = lo->cor[idx_d];
    }

    // --->
    //      <---
    // The current lo is next to, and not overlapping, current hi
    else
    {
      hi->adj[!d] = lo;
      hi = hi->cor[idx_c];
    }
  }
}
void link_adjacent_blocks(Synmap * syn){
  for (size_t genid = 0; genid <= 1; genid++){
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
  Node * root;
  long qdiff, tdiff;
  char this_strand, older_strand;
  size_t setid = 0;
  for (size_t i = 0; i < SG(syn,0)->size; i++) {
    // Initialize the first block in the scaffold
    blk = SGC(syn, 0, i)->head[0]; 
    blk->setid = ++setid;
    blk->over->setid = blk->setid;
    node = init_node(blk);
    root = node;
    for (blk = blk->cor[1]; blk != NULL; blk = blk->cor[1]){
      this_strand = blk->over->strand;
      while(true){
        // qdiff and tdiff describe the adjacency of blocks relative to the
        // query are target contigs, respectively. Cases:
        // ---
        // diff <= -2 : blocks are not adjacent
        // diff == -1 : blocks are adjacent on reverse strand
        // diff ==  0 : blocks overlap
        // diff ==  1 : blocks are adjacent
        // diff >=  2 : blocks are not adjacent
        qdiff = (long) blk->grpid       - (long) node->blk->grpid;
        tdiff = (long) blk->over->grpid - (long) node->blk->over->grpid;

        older_strand = node->blk->over->strand;

        // If adjacent
        if((qdiff == 1) &&
           (this_strand == older_strand) &&
           ((tdiff ==  1 && this_strand == '+') ||
            (tdiff == -1 && this_strand == '-')))
        {
          blk->setid = node->blk->setid;
          blk->over->setid = blk->setid;

          blk->cnr[0] = node->blk;
          blk->over->cnr[0] = node->blk->over;

          node->blk->cnr[1] = blk;
          node->blk->over->cnr[1] = blk->over;

          node->blk = blk;

          // TODO: in strange cases, one node might be adjacent to multiple
          // nodes. By placing a break here, I just take the first. I need
          // explicit handling for this case.
          break;
        }
        // If at bottom
        else if(node->down == NULL){
          blk->setid = ++setid;
          blk->over->setid = blk->setid;
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
    free_node(root);
  }
}

// TODO extend validatation to corners and the other new stuff
void validate_synmap(Synmap * syn){
    size_t gid, cid;
    Contig * con;
    Block  * blk;
    assert(syn->size == 2);
    for(gid = 0; gid < syn->size; gid++){
        for(cid = 0; cid < SG(syn, gid)->size; cid++){
            con = SGC(syn, gid, cid);
            blk = con->head[0];
            for(; blk != NULL; blk = blk->cor[1]){
                // assert(blk->pos[1] < con->length);
                if(!(blk->pos[1] < con->length)){
                    fprintf(stderr,
                            "WARNING: stop greater than contig length: %zu vs %zu\n",
                            blk->pos[1], con->length);
                }
                assert(blk->setid == blk->over->setid);
                assert(blk->score == blk->over->score);
                assert(blk->setid != 0);
                assert(blk->grpid != 0);
                if(blk->cnr[1] != NULL){
                    assert(blk->grpid != blk->cnr[1]->grpid);
                    assert(blk->setid == blk->cnr[1]->setid);
                    assert(blk->cnr[1]->over->cnr[0] != NULL);
                    assert(blk->cnr[1]->over->cnr[0]->over == blk);
                }
            }
        }
    }
}
