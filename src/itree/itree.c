#include "itree.h"

/* local function prototypes */
IntervalTree * build_tree_r(IA * intervals, IntervalTree * parent, Orientation orientation);
void print_IntervalTree_verbosity_1(IntervalTree * n, int depth, char pos);
void print_IntervalTree_verbosity_2(IntervalTree * n, int depth, char pos);
void print_IntervalTree_verbosity_3(IntervalTree * n, int depth, char pos);
void print_IntervalTree_r(IntervalTree * n, int depth, char pos, int verbosity);
long get_center(IA *);


IntervalTree * init_IntervalTree(){
    IntervalTree * tree = (IntervalTree *)malloc(sizeof(IntervalTree));
    tree->center   = 0;
    tree->by_start = NULL;
    tree->by_stop  = NULL;
    tree->children[0] = NULL;
    tree->children[1] = NULL;
    tree->parent = NULL;
    tree->orientation = O_UNSET;
    return(tree);
}

void free_IntervalTree(IntervalTree * tree){
    if(LEFT(tree))
        free_IntervalTree(LEFT(tree));
    if(RIGHT(tree))
        free_IntervalTree(RIGHT(tree));
    if(tree->by_start)
        free_IA(tree->by_start);
    if(tree->by_stop)
        free_IA(tree->by_stop);
    if(tree != NULL)
        free(tree);
}

IntervalTree * build_tree(IA * intervals){
    return build_tree_r(intervals, NULL, O_ROOT);
}

IntervalTree * build_tree_r(IA * intervals, IntervalTree * parent, Orientation orientation){
    /* initialize returned product */
    IntervalTree * tree = init_IntervalTree();
    tree->parent = parent;
    tree->orientation = orientation;
    tree->center = get_center(intervals);

    /* array to store position of center point relative to each interval
     * in intervals
     * lo = 0 -> center is lower than interval
     * in = 1 -> center is inside interval
     * hi = 2 -> center is higher than interval
     */
    Pos * pos = (Pos *)malloc(intervals->size * sizeof(Pos));

    /* count of intervals in each relative position */
    int npos[] = {0, 0, 0};

    /* iterate over intervals classifying and counting each */
    for(size_t i = 0; i < intervals->size; i++){
        pos[i] = point_overlap(tree->center, &intervals->v[i]);
        npos[pos[i]]++;
    }

    /* initialise interval arrays */
    IA * right     = init_set_IA(npos[lo]);
    IA * left      = init_set_IA(npos[hi]);
    tree->by_start = init_set_IA(npos[in]);
    tree->by_stop  = init_set_IA(npos[in]);

    /* track index of interval to add */
    size_t lo_idx = 0;
    size_t in_idx = 0;
    size_t hi_idx = 0;

    /* assign intervals to appropriate partitions */
    for(size_t i = 0; i < intervals->size; i++){
        switch(pos[i]){
            case lo:
                right->v[lo_idx] = intervals->v[i];
                lo_idx++;
                break;
            case in:
                tree->by_start->v[in_idx] = intervals->v[i];
                tree->by_stop->v[in_idx]  = intervals->v[i];
                in_idx++;
                break;
            case hi:
                left->v[hi_idx] = intervals->v[i];
                hi_idx++;
                break;
            default:
                fprintf(stderr, "ERROR: something really weird happened ... sorry\n");
                exit(EXIT_FAILURE);
        }
    }

    if(npos[lo] > 0){
        RIGHT(tree) = build_tree_r(right, tree, O_RIGHT);
    } else {
        free_IA(right);
    }

    if(npos[hi] > 0){
        LEFT(tree) = build_tree_r(left, tree, O_LEFT);
    } else {
        free_IA(left);
    }
    
    if(npos[in] > 0){
        qsort(&tree->by_start->v[0], npos[in], sizeof(Interval), cmp_start);
        qsort(&tree->by_stop->v[0],  npos[in], sizeof(Interval), cmp_stop);
    }

    free(pos);
    free_IA(intervals);

    return(tree);
}



/*
 * Select a point at the center of the middle interval.
 * This guarantees at least one interval overlaps each node.
 * If the intervals are sorted, it also favors (but doesn't guarantee) a
 * balanced tree.
 */
long get_center(IA * intr){
    // get the central index
    long i = intr->size / 2;
    // get the center point on this index
    long x = (intr->v[i].stop - intr->v[i].start) / 2 + intr->v[i].start;
    return x;
}

/* write tree and center */
void print_IntervalTree_verbosity_1(IntervalTree * n, int depth, char pos){
    printf("%*d - %c%zu\n", depth * 2, depth, pos, n->center);
}

/* write tree, center, and start-sorted */
void print_IntervalTree_verbosity_2(IntervalTree * n, int depth, char pos){
    printf("%*d   %*s\t%c%zu:",
           depth * 2, depth, 
           10 - depth * 2, "", pos, n->center);
    for(size_t i = 0; i < n->by_start->size; i++){
        printf("(%zu,%zu) ",
               n->by_start->v[i].start,
               n->by_start->v[i].stop);
    }
    printf("\n");
}

/* write start- and stop-sorted vectors for each node */
void print_IntervalTree_verbosity_3(IntervalTree * n, int depth, char pos){
    print_IntervalTree_verbosity_1(n, depth, pos);
    for(size_t i = 0; i < n->by_start->size; i++){
        printf("\t\t(%zu,%zu) ",
               n->by_start->v[i].start,
               n->by_start->v[i].stop);
        printf("(%zu,%zu)\n",
               n->by_stop->v[i].start,
               n->by_stop->v[i].stop);
    }
}

/* local print function */
void print_IntervalTree_r(IntervalTree * n, int depth, char pos, int verbosity){
    switch(verbosity){
        case 1:
            print_IntervalTree_verbosity_1(n, depth, pos); break;
        case 2:
            print_IntervalTree_verbosity_2(n, depth, pos); break;
        case 3:
            print_IntervalTree_verbosity_3(n, depth, pos); break;
        default:
            fprintf(stderr, "verbosity must be 1, 2, or 3\n");
            exit(EXIT_FAILURE);
    }
    depth++;
    if(LEFT(n) != NULL){
        print_IntervalTree_r(LEFT(n), depth, 'l', verbosity);
    }
    if(RIGHT(n) != NULL){
        print_IntervalTree_r(RIGHT(n), depth, 'r', verbosity);
    }
}

/* public wrapper for real print function */
void print_IntervalTree(IntervalTree * n, int verbosity){
    print_IntervalTree_r(n, 0, 'c', verbosity);
}
