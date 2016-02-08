#include <stdlib.h>

#include "interval-tree.h"
#include "ia.h"

struct IntervalTree * init_interval_tree(){
    struct IntervalTree * tree = (struct IntervalTree *)malloc(sizeof(struct IntervalTree));
    tree->center   = 0;
    tree->by_start = NULL;
    tree->by_stop  = NULL;
    tree->l_child  = NULL;
    tree->r_child  = NULL;
    return(tree);
}
void free_interval_tree(struct IntervalTree * tree){
    if(tree->l_child)
        free_interval_tree(tree->l_child);
    if(tree->r_child)
        free_interval_tree(tree->r_child);
    if(tree->by_start)
        free_ia(tree->by_start);
    if(tree->by_stop)
        free_ia(tree->by_stop);
    if(tree)
        free(tree);
}
