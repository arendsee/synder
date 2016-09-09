#include "itree_node.h"

IntervalTreeNode::IntervalTreeNode(
    std::vector<Interval*> intervals,
    IntervalTreeNode * new_parent,
    Orientation new_orientation
) :
    parent(new_parent), orientation(new_orientation)
{

    center = get_center(intervals);

    std::vector<Interval*> left;
    std::vector<Interval*> right;

    /* iterate over intervals classifying and counting each */
    Interval * v;
    for(size_t i = 0; i < intervals.size(); i++){
        v = intervals[i];
        switch(v->overlap(center)){
            case lo:
                right.push_back(v);
                break;
            case hi:
                left.push_back(v);
                break;
            case in:
                by_stop.push_back(v);
                by_start.push_back(v);
                break;
            default:
                break;
        }
    }

    if(!by_start.empty())
        std::sort(by_start.begin(), by_start.end(), Interval::cmp_start);

    if(!by_stop.empty())
        std::sort(by_stop.begin(), by_stop.end(), Interval::cmp_stop);

    if (left.empty()){
        children[0] = NULL;
    } else {
        children[0] = new IntervalTreeNode(left, this, O_LEFT);
    }

    if (right.empty()){
        children[1] = NULL;
    } else {
        children[1] = new IntervalTreeNode(right, this, O_RIGHT);
    }
}


IntervalTreeNode::~IntervalTreeNode(){
    if (children[0] != NULL)
        delete children[0];

    if (children[1] != NULL)
        delete children[1];
}


/*
 * Select a point at the center of the middle interval.
 * This guarantees at least one interval overlaps each node.
 * If the intervals are sorted, it also favors (but doesn't guarantee) a
 * balanced tree.
 */
long IntervalTreeNode::get_center(std::vector<Interval*> intr){
    // get the central index
    long i = intr.size() / 2;
    // get the center point on this index
    long x = (intr[i]->stop - intr[i]->start) / 2 + intr[i]->start;
    return x;
}

/* write tree and center */
void IntervalTreeNode::print_verbosity_1(IntervalTreeNode * n, int depth, char pos){
    printf("%*d - %c%zu\n", depth * 2, depth, pos, n->center);
}

/* write tree, center, and start-sorted */
void IntervalTreeNode::print_verbosity_2(IntervalTreeNode * n, int depth, char pos){
    printf("%*d   %*s\t%c%zu:",
           depth * 2, depth, 
           10 - depth * 2, "", pos, n->center);
    for(size_t i = 0; i < n->by_start.size(); i++){
        printf("(%zu,%zu) ",
               n->by_start[i]->start,
               n->by_start[i]->stop);
    }
    printf("\n");
}

/* write start- and stop-sorted vectors for each node */
void IntervalTreeNode::print_verbosity_3(IntervalTreeNode * n, int depth, char pos){
    print_verbosity_1(n, depth, pos);
    for(size_t i = 0; i < n->by_start.size(); i++){
        printf("\t\t(%zu,%zu) ",
               n->by_start[i]->start,
               n->by_start[i]->stop);
        printf("(%zu,%zu)\n",
               n->by_stop[i]->start,
               n->by_stop[i]->stop);
    }
}

/* local print function */
void IntervalTreeNode::print(IntervalTreeNode * n, int depth, char pos, int verbosity){
    switch(verbosity){
        case 1:
            print_verbosity_1(n, depth, pos); break;
        case 2:
            print_verbosity_2(n, depth, pos); break;
        case 3:
            print_verbosity_3(n, depth, pos); break;
        default:
            fprintf(stderr, "verbosity must be 1, 2, or 3\n");
            exit(EXIT_FAILURE);
    }
    depth++;
    if(LEFT(n) != NULL){
        print(LEFT(n), depth, 'l', verbosity);
    }
    if(RIGHT(n) != NULL){
        print(RIGHT(n), depth, 'r', verbosity);
    }
}

/* public wrapper for real print function */
void IntervalTreeNode::print(int verbosity){
    print(this, 0, 'c', verbosity);
}
