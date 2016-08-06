void print_ContiguousMap(ContiguousMap * cmap){
  for(int i = 0; i < cmap->size; i++){
    printf("\n--- NODE %i ---\n", i);
    print_ContiguousNode(cmap->map[i]);
  }
}

void print_ContiguousNode(ContiguousNode * cnode){
  printf("feature:\n");
  print_Block(cnode->feature);
  printf("match:\n");
  print_Block(cnode->match);
  printf("prev: %p\n", cnode->prev);
  printf("next: %p\n", cnode->next);
  printf("flag=%i;qblkid=%lu;setid=%lu\n",
         cnode->flag, cnode->qblkid, cnode->setid);
}


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


// ---- A local utility data structure used to build contiguous sets ----
// INVARIANT: root should never change after initialization
typedef struct NodeList{
  struct NodeList * down;
  struct NodeList * up;
  ContiguousNode * cnode; 
} NodeList;

NodeList * init_NodeList(ContiguousNode * cnode){
  NodeList * node = (NodeList *)malloc(sizeof(NodeList));
  node->cnode = cnode;
  node->down = NULL;
  node->up = NULL;
  return(node);
}

void init_down_NodeList(NodeList * node, ContiguousNode * cnode){
  if(node->cnode == NULL){
    node->cnode = cnode; 
  } else {
    node->down = init_NodeList(cnode);
    node->down->up = node;
  }
}

NodeList * remove_node_NodeList(NodeList * node){
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

void free_NodeList(NodeList * node){
  if(node->down != NULL){
    free_NodeList(node->down);
  }
  free(node);
}

void print_NodeList(NodeList * node, int level){
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
