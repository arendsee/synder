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

void link_block_corners(Synmap* syn)
{
    Contig* con;
    size_t N;
    for (size_t g = 0; g <= 1; g++) {
        for (size_t c = 0; c < SG(syn, g)->size; c++) {
            con = SGC(syn, g, c);
            N = con->size;

            Block ** blocks = (Block **)malloc(N * sizeof(Block *));
            for(size_t i = 0; i < N; i++){
                blocks[i] = &con->block[i];
            }
            
            // sort by stop
            sort_blocks(blocks, N, true);
            for (size_t i = 0; i < N; i++) {
                blocks[i]->cor[PREV_STOP] = (i == 0)     ? NULL : blocks[i - 1];
                blocks[i]->cor[NEXT_STOP] = (i == N - 1) ? NULL : blocks[i + 1];
            }
            // sort by start
            sort_blocks(blocks, N, false);
            for (size_t i = 0; i < N; i++) {
                blocks[i]->cor[PREV_START] = (i == 0)     ? NULL : blocks[i - 1];
                blocks[i]->cor[NEXT_START] = (i == N - 1) ? NULL : blocks[i + 1];
            }

            free(blocks);
        }
    }
}

void set_contig_corners(Synmap * syn){
    Contig * con;
    for(size_t g = 0; g < 2; g++){
        for(size_t c = 0; c < SG(syn, g)->size; c++){
            con = SGC(syn, g, c); 
            con->cor[1] = &con->block[0];
            con->cor[0] = &con->block[0];
            con->cor[2] = &con->block[con->size-1];
            con->cor[3] = &con->block[con->size-1];
            while(con->cor[1]->cor[PREV_STOP] != NULL) {
                con->cor[1] = con->cor[1]->cor[PREV_STOP];
            }
            while(con->cor[0]->cor[PREV_START] != NULL) {
                con->cor[0] = con->cor[0]->cor[PREV_START];
            }
            while(con->cor[2]->cor[NEXT_START] != NULL) {
                con->cor[2] = con->cor[2]->cor[NEXT_STOP];
            }
            while(con->cor[3]->cor[NEXT_STOP] != NULL) {
                con->cor[3] = con->cor[3]->cor[NEXT_STOP];
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
      blk = SGC(syn,g,i)->cor[0];
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

  if(con->cor[d] == NULL || con->cor[!d] == NULL){
    fprintf(stderr, "Contig head must be set before link_adjacent_blocks is called\n");
    exit(EXIT_FAILURE);
  }

  // Transformed indices for Block->cor
  // ----------------------------------
  // a   b          c   d
  // <--->   <--->  <--->
  //   a - previous element by start
  //   b - previous element by stop
  //   c - next element by start
  //   d - next element by stop
  int idx_a = (!d * 2) + !d;  // - 0
  // int idx_b = ( d * 2) + !d;  // - 2
  int idx_c = (!d * 2) +  d;  // - 1
  int idx_d = ( d * 2) +  d;  // - 3

  Block *lo, *hi;

  lo = con->cor[idx_c]; // first element by stop
  hi = con->cor[idx_a]; // first element by start
  
  while(hi != NULL){

    //       ---> 
    // <---
    // OR
    // --->
    //   <---
    // This should should occur only at the beginning
    if(REL_LE(hi->pos[!d], lo->pos[d], d))
    {
      hi->adj[!d] = NULL;
      hi = hi->cor[idx_c];
    }

    //  lo     next
    // ---->  ---->
    //               <---
    // If next is closer, and not overlapping the hi, increment lo
    // You increment lo until it is adjacent to the current hi
    else if(REL_LT(lo->cor[idx_d]->pos[d], hi->pos[!d], d))
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
  node->down  = NULL;
  node->blk   = blk;
  return(node);
}
void remove_node(Node * node)
{
  if(node->down != NULL){
    Node * tmp = node->down;
    memcpy(node, node->down, sizeof(Node));
    free(tmp);
  } else {
    node->blk = NULL;
  }
}
void free_node(Node * node)
{
  if(node->down != NULL)
    free_node(node->down);
  free(node);
}

void link_contiguous_blocks(Synmap * syn, long k)
{
  Block * blk;
  Node * node;
  Node * root;
  long qdiff, tdiff;
  char this_strand, older_strand;
  size_t setid = 0;
  bool same_target;
  for (size_t i = 0; i < SG(syn,0)->size; i++) {
    // Initialize the first block in the scaffold
    blk = SGC(syn, 0, i)->cor[0]; 
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
        same_target = blk->over->parent == node->blk->over->parent;

        // If contiguous
        if(
              // query blocks are contiguous
              (
                  // query blocks are adjacent
                  (qdiff == 1) ||
                  // fewer than cutoff blocks with different targets separate them
                  (same_target && (qdiff-1) < k && qdiff != 0)
              ) &&
              // targets are contiguous
              (
                  // target intervals are on the same strand
                  this_strand == older_strand &&
                  // target intervals are not overlapping 
                  tdiff != 0 &&
                  (
                      // targets are immediately adjacent
                      (
                          (tdiff ==  1 && this_strand == '+')  ||
                          (tdiff == -1 && this_strand == '-')
                      ) ||
                      // targets are close enough
                      (
                          // target stuff is in the right direction
                          (
                              (tdiff > 0 && this_strand == '+') ||
                              (tdiff < 0 && this_strand == '-') 
                          )
                          && // score is acceptable
                          ((abs(tdiff) + qdiff - 2) <= k)
                      )
                  )
              )
          )
        {
          blk->setid = node->blk->setid;
          blk->over->setid = blk->setid;

          blk->cnr[0] = node->blk;
          blk->over->cnr[0] = node->blk->over;

          node->blk->cnr[1] = blk;
          node->blk->over->cnr[1] = blk->over;

          node->blk = blk;

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
        else if((qdiff - 1) > k){
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
    size_t nblks;
    Contig * con;
    Block  * blk;
    assert(syn->size == 2);
    for(gid = 0; gid < syn->size; gid++){
        for(cid = 0; cid < SG(syn, gid)->size; cid++){
            con = SGC(syn, gid, cid);
            blk = con->cor[0];
            nblks = 0;
            for(; blk != NULL; blk = blk->cor[1]){
                if(!(blk->pos[1] < con->length)){
                    fprintf(stderr,
                            "WARNING: stop greater than contig length: %zu vs %zu\n",
                            blk->pos[1], con->length);
                }
                nblks++;
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
                if(blk->cor[NEXT_START] != NULL){
                    assert(blk->pos[0] <= blk->cor[NEXT_START]->pos[0]);
                    assert(blk->cor[NEXT_START]->cor[PREV_START] != NULL);
                    assert(blk == blk->cor[NEXT_START]->cor[PREV_START]);
                }
                if(blk->cor[NEXT_STOP] != NULL){
                    assert(blk->pos[1] <= blk->cor[NEXT_STOP]->pos[1]);
                    assert(blk->cor[NEXT_STOP]->cor[PREV_STOP] != NULL);
                    assert(blk == blk->cor[NEXT_STOP]->cor[PREV_STOP]);
                }
                if(blk->cor[PREV_START] != NULL){
                    assert(blk->cor[PREV_START]->cor[NEXT_START] != NULL);
                }
                if(blk->cor[PREV_STOP] != NULL){
                    assert(blk->cor[PREV_STOP]->cor[NEXT_STOP] != NULL);
                }
            }
            assert(nblks <= con->size);
        }
    }
}
