#include "search.hpp"

long count_interval_overlaps_r(Interval *, IntervalTree *, long);
long count_point_overlaps_r(long, IntervalTree *, long);

IntervalResult * get_interval_overlaps_r(Interval *, IntervalTree *, IntervalResult *);
void get_point_overlaps_r(long, IntervalTree *, IntervalResult *);

long count_interval_overlaps(Interval * inv, IntervalTree * tree){
    return count_interval_overlaps_r(inv, tree, 0);
}

long count_point_overlaps(long point, IntervalTree * tree){
    return count_point_overlaps_r(point, tree, 0);
}

IntervalResult * init_IntervalResult(){
    IntervalResult * res = (IntervalResult *)malloc(sizeof(IntervalResult));
    res->iv = NULL;
    res->inbetween = false;
    res->leftmost = false;
    res->rightmost = false;
    return(res);
}

void print_IntervalResult(IntervalResult * res){
    print_IV(res->iv);
    printf("inbetween=%i leftmost=%i rightmost=%i\n", res->inbetween, res->leftmost, res->rightmost);
}

void free_IntervalResult(IntervalResult * res){
    if(res != NULL){
        if(res->iv != NULL)
            free_IV(res->iv);
        free(res);
    }
}

IntervalResult * get_point_overlaps(long point, IntervalTree * tree){
    IntervalResult * res = init_IntervalResult();
    res->iv = init_IV(IV_INITIAL_SIZE); 
    get_point_overlaps_r(point, tree, res);
    return res;
}

IntervalResult * get_interval_overlaps(Interval * inv, IntervalTree * tree){
    IntervalResult * res = init_IntervalResult();
    res->iv = init_IV(IV_INITIAL_SIZE); 
    get_interval_overlaps_r(inv, tree, res);
    return res;
}

long count_point_overlaps_r(long point, IntervalTree * tree, long count){
    if(point >= tree->center) {
        for(long long i = T_SIZE(tree) - 1; i >= 0 ; i--){
            if(point <= T_STOP_STOP(tree, i)){
                count++;
            } else {
                break;
            }
        }
        if(RIGHT(tree) != NULL)
            return count_point_overlaps_r(point, RIGHT(tree), count);
    }
    else {
        for(size_t i = 0; i < T_SIZE(tree); i++){
            if(point >= T_START_START(tree, i)){
                count++;
            } else {
                break;
            }
        }
        if(LEFT(tree) != NULL)
            return count_point_overlaps_r(point, LEFT(tree), count);
    }
    return count;
}

long count_interval_overlaps_r(Interval * inv, IntervalTree * tree, long count){
    if(tree == NULL)
        return count;
    switch(point_overlap(tree->center, inv)){
    case lo:
        for(long long i = T_SIZE(tree) - 1; i >= 0 ; i--){
            if(inv->start <= T_STOP_STOP(tree, i)){
                count++;
            } else {
                break;
            }
        }
        return count_interval_overlaps_r(inv, RIGHT(tree), count);
    case hi:
        for(size_t i = 0; i < T_SIZE(tree); i++){
            if(inv->stop >= T_START_START(tree, i)){
                count++;
            } else {
                break;
            }
        }
        return count_interval_overlaps_r(inv, LEFT(tree), count);
    default:
        count += T_SIZE(tree);
        return count_interval_overlaps_r(inv, RIGHT(tree),
               count_interval_overlaps_r(inv, LEFT(tree), count));
    }
}

void set_nearest_opposing_interval(IntervalTree * node, IntervalResult * results, Pos pos){
    
    Orientation orientation = node->orientation;

    if((pos == hi && orientation == O_LEFT) ||
       (pos == lo && orientation == O_RIGHT))
    {
        node = node->parent;    
    }
    else {
        while(node->orientation != O_ROOT && node->orientation == orientation){
            node = node->parent;
        }
        node = node->parent;
    }

    if(node == NULL){
      if(pos == lo){
        results->leftmost = true;
      } else{
        results->rightmost = true;
      }
    } else {
        results->inbetween = true;
        if(orientation == O_RIGHT){
          add_IV(results->iv, node->by_start->v[0]);
        } else if(orientation == O_LEFT) {
          add_IV(results->iv, node->by_stop->v[node->by_stop->size - 1]);
        }
    }
}

void get_point_overlaps_r(long point, IntervalTree * tree, IntervalResult * results){
    if(point >= tree->center) {
        for(long long i = T_SIZE(tree) - 1; i >= 0 ; i--){
            if(point <= T_STOP_STOP(tree, i)){
                add_IV(results->iv, T_STOP(tree, i));
            } else {
                break;
            }
        }
        if(RIGHT(tree) != NULL){
            get_point_overlaps_r(point, RIGHT(tree), results);
        } 
        else if(R_SIZE(results) != 0) {
            add_IV(results->iv, LAST_STOP(tree));
            set_nearest_opposing_interval(tree, results, hi);
        }
    }
    else {
        for(size_t i = 0; i < T_SIZE(tree); i++){
            if(point >= T_START_START(tree, i)){
                add_IV(results->iv, T_START(tree, i));
            } else {
                break;
            }
        }
        if(LEFT(tree) != NULL){
            get_point_overlaps_r(point, LEFT(tree), results);
        }
        else if(R_SIZE(results) != 0) {
            add_IV(results->iv, FIRST_START(tree));
            set_nearest_opposing_interval(tree, results, lo);
        }
    }
}

IntervalResult * get_interval_overlaps_r(Interval * inv, IntervalTree * tree, IntervalResult * results){
    if(tree == NULL)
        return results;
    switch(point_overlap(tree->center, inv)){
    case lo: // center lower than interval start
        for(long long i = T_SIZE(tree) - 1; i >= 0 ; i--){
            if(inv->start <= T_STOP_STOP(tree, i)){
                add_IV(results->iv, T_STOP(tree, i));
            } else {
                break;
            }
        }
        // Reach a leaf and still have no overlaps
        if(RIGHT(tree) == NULL && results->iv->size == 0){
            add_IV(results->iv, LAST_STOP(tree));
            // get nearest interval on the other side
            set_nearest_opposing_interval(tree, results, hi);
        }
        return get_interval_overlaps_r(inv, RIGHT(tree), results);
    case hi:
        for(size_t i = 0; i < T_SIZE(tree); i++){
            if(inv->stop >= T_START_START(tree, i)){
                add_IV(results->iv, T_START(tree, i));
            } else {
                break;
            }
        }
        // Reach a leaf and still have no overlaps
        if(LEFT(tree) == NULL && results->iv->size == 0){
            add_IV(results->iv, FIRST_START(tree));
            // get nearest interval on the other side
            set_nearest_opposing_interval(tree, results, lo);
        }
        return get_interval_overlaps_r(inv, LEFT(tree), results);
    default: // in
        join_IV(results->iv, tree->by_start);
        return get_interval_overlaps_r(inv, RIGHT(tree), 
               get_interval_overlaps_r(inv, LEFT(tree), results));
    };
}
