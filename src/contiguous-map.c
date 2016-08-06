#include "contiguous-map.h"

// ---- A local utility structure used to build contiguous sets ----
// INVARIANT: root should never change after initialization
typedef struct NodeList{
  struct NodeList * down;
  struct NodeList * up;
  ContiguousNode * cnode; 
} NodeList;
NodeList * init_NodeList(ContiguousNode * cnode);
void print_NodeList(NodeList * node, int level);
void free_NodeList(NodeList * node);
NodeList * remove_node_NodeList(NodeList * node);
void init_down_NodeList(NodeList * node, ContiguousNode * cnode);


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



ContiguousNode * init_ContiguousNode(Synmap * syn, size_t conid, size_t blkid)
{
  ContiguousNode *cnode = (ContiguousNode *) malloc(sizeof(ContiguousNode));
  cnode->feature = syn->genome[0]->contig[conid]->block[blkid];
  cnode->next = NULL;
  cnode->prev = NULL;
  cnode->match = QT_SGCB(syn, cnode->feature);
  cnode->qblkid = blkid;
  return cnode;  
}

ContiguousMap *init_ContiguousMap(size_t size)
{
  ContiguousMap *cmap = (ContiguousMap *) malloc(sizeof(ContiguousMap));
  cmap->size = size;
  cmap->map = (ContiguousNode **) malloc(size * sizeof(ContiguousNode *));
  return cmap;
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

void print_ContiguousMap(ContiguousMap * cmap)
{
  for(int i = 0; i < cmap->size; i++){
    printf("\n--- NODE %i ---\n", i);
    print_ContiguousNode(cmap->map[i]);
  }
}

void print_ContiguousNode(ContiguousNode * cnode)
{
  printf("feature:\n");
  print_Block(cnode->feature);
  printf("match:\n");
  print_Block(cnode->match);
  printf("prev: %p\n", cnode->prev);
  printf("next: %p\n", cnode->next);
  printf("qblkid=%lu;setid=%lu\n", cnode->qblkid, cnode->setid);
}

NodeList * init_NodeList(ContiguousNode * cnode)
{
  NodeList * node = (NodeList *)malloc(sizeof(NodeList));
  node->cnode = cnode;
  node->down = NULL;
  node->up = NULL;
  return(node);
}

void init_down_NodeList(NodeList * node, ContiguousNode * cnode)
{
  if(node->cnode == NULL){
    node->cnode = cnode; 
  } else {
    node->down = init_NodeList(cnode);
    node->down->up = node;
  }
}

NodeList * remove_node_NodeList(NodeList * node)
{
  NodeList * tmp;
  if(node->up == NULL){ // is root
    if(node->down == NULL){
      node->cnode = NULL;
    } else {
      tmp = node->down;
      if(node->down->down != NULL){
        node->down->down->up = node;
      }
      node->cnode = node->down->cnode;
      node->down = node->down->down;
      free(tmp);
    }
  } else {
    tmp = node;
    node->up->down = node->down;
    if(node->down != NULL){
      node->down->up = node->up;
    }
    node = node->up;
    free(tmp);
  }
  return node;
}

void free_NodeList(NodeList * node)
{
  if(node->down != NULL){
    free_NodeList(node->down);
  }
  free(node);
}

void print_NodeList(NodeList * node, int level)
{
  if(node->cnode != NULL){
    printf(
      "%i (%u, %u)\n",
      level,
      node->cnode->feature->start,
      node->cnode->feature->stop
    );
  }
  if(node->down != NULL){
    print_NodeList(node->down, level++);
  }
}
