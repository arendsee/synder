#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

#include "search.h"
#include "iv.h"

uint count_interval_overlaps_r(Interval *, IntervalTree *, uint);
uint count_point_overlaps_r(uint, IntervalTree *, uint);

IntervalResult * get_interval_overlaps_r(Interval *, IntervalTree *, IntervalResult *);
void get_point_overlaps_r(uint, IntervalTree *, IntervalResult *);

uint count_interval_overlaps(Interval * inv, IntervalTree * tree){
    return count_interval_overlaps_r(inv, tree, 0);
}

uint count_point_overlaps(uint point, IntervalTree * tree){
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
    print_iv(res->iv);
    printf("inbetween=%i leftmost=%i rightmost=%i\n", res->inbetween, res->leftmost, res->rightmost);
}

void free_IntervalResult(IntervalResult * res){
    if(res != NULL){
        if(res->iv != NULL)
            iv_free(res->iv);
        free(res);
    }
}

IntervalResult * get_point_overlaps(uint point, IntervalTree * tree){
    IntervalResult * res = init_IntervalResult();
    res->iv = iv_init(IV_INITIAL_SIZE); 
    get_point_overlaps_r(point, tree, res);
    return res;
}

IntervalResult * get_interval_overlaps(Interval * inv, IntervalTree * tree){
    IntervalResult * res = init_IntervalResult(IV_INITIAL_SIZE);
    res->iv = iv_init(IV_INITIAL_SIZE); 
    get_interval_overlaps_r(inv, tree, res);
    return res;
}

uint count_point_overlaps_r(uint point, IntervalTree * tree, uint count){
    if(point >= tree->center) {
        for(int i = T_SIZE(tree) - 1; i >= 0 ; i--){
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
        for(int i = 0; i < T_SIZE(tree); i++){
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

uint count_interval_overlaps_r(Interval * inv, IntervalTree * tree, uint count){
    if(tree == NULL)
        return count;
    switch(point_overlap(tree->center, inv)){
    case lo:
        for(int i = T_SIZE(tree) - 1; i >= 0 ; i--){
            if(inv->start <= T_STOP_STOP(tree, i)){
                count++;
            } else {
                break;
            }
        }
        return count_interval_overlaps_r(inv, RIGHT(tree), count);
    case hi:
        for(int i = 0; i < T_SIZE(tree); i++){
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
    /*
     * Two nodes are CIS if they have the same orientation, e.g. both are left nodes.
     *
     *                   a
     *                /     \
     *               b       e
     *              / \     / \
     *             c   d   f   g
     *
     *  c is CIS to b and f
     *  d is CIS to e and g
    */
    #define SET_INITIAL_ORIENTATION \
      Orientation orientation = node->orientation;
    #define ASCEND            node = node->parent;
    #define IS_ROOT           node->orientation == O_ROOT
    #define IS_CIS            node->orientation == orientation
    #define SET_EXTREME \
      if(pos == lo){ results->leftmost = true; } \
      else { results->rightmost = true; }
    #define SET_NEAREST_NODE \
      if(orientation == O_RIGHT) \
      { iv_add(results->iv, node->by_start->v[0]); } \
      else if(orientation == O_LEFT) \
      { iv_add(results->iv, node->by_stop->v[node->by_stop->size - 1]); }

    SET_INITIAL_ORIENTATION 

    if((pos == hi && orientation == O_LEFT) ||
       (pos == lo && orientation == O_RIGHT))
    {
        ASCEND    
    }
    else {
        while(!IS_ROOT && IS_CIS){ ASCEND }
        ASCEND
    }

    if(node == NULL){
        SET_EXTREME
    } else {
        SET_NEAREST_NODE 
    }

    #undef ASCEND
    #undef IS_ROOT
    #undef IS_CIS
    #undef SET_EXTREME
}

void get_point_overlaps_r(uint point, IntervalTree * tree, IntervalResult * results){
    if(point >= tree->center) {
        for(int i = T_SIZE(tree) - 1; i >= 0 ; i--){
            if(point <= T_STOP_STOP(tree, i)){
                iv_add(results->iv, T_STOP(tree, i));
            } else {
                break;
            }
        }
        if(RIGHT(tree) != NULL){
            get_point_overlaps_r(point, RIGHT(tree), results);
        } 
        else if(R_SIZE(results) != 0) {
            results->inbetween = true;
            iv_add(results->iv, LAST_STOP(tree));
            set_nearest_opposing_interval(tree, results, hi);
        }
    }
    else {
        for(int i = 0; i < T_SIZE(tree); i++){
            if(point >= T_START_START(tree, i)){
                iv_add(results->iv, T_START(tree, i));
            } else {
                break;
            }
        }
        if(LEFT(tree) != NULL){
            get_point_overlaps_r(point, LEFT(tree), results);
        }
        else if(R_SIZE(results) != 0) {
            results->inbetween = true;
            iv_add(results->iv, FIRST_START(tree));
            set_nearest_opposing_interval(tree, results, lo);
        }
    }
}

IntervalResult * get_interval_overlaps_r(Interval * inv, IntervalTree * tree, IntervalResult * results){
    if(tree == NULL)
        return results;
    switch(point_overlap(tree->center, inv)){
    case lo: // center lower than interval start
        for(int i = T_SIZE(tree) - 1; i >= 0 ; i--){
            if(inv->start <= T_STOP_STOP(tree, i)){
                iv_add(results->iv, T_STOP(tree, i));
            } else {
                break;
            }
        }
        // Reach a leaf and still have no overlaps
        if(RIGHT(tree) == NULL && results->iv->size == 0){
            iv_add(results->iv, LAST_STOP(tree));
            // get nearest interval on the other side
            set_nearest_opposing_interval(tree, results, hi);
        }
        return get_interval_overlaps_r(inv, RIGHT(tree), results);
    case hi:
        for(int i = 0; i < T_SIZE(tree); i++){
            if(inv->stop >= T_START_START(tree, i)){
                iv_add(results->iv, T_START(tree, i));
            } else {
                break;
            }
        }
        // Reach a leaf and still have no overlaps
        if(LEFT(tree) == NULL && results->iv->size == 0){
            iv_add(results->iv, FIRST_START(tree));
            // get nearest interval on the other side
            set_nearest_opposing_interval(tree, results, lo);
        }
        return get_interval_overlaps_r(inv, LEFT(tree), results);
    default: // in
        iv_join(results->iv, tree->by_start);
        return get_interval_overlaps_r(inv, RIGHT(tree), 
               get_interval_overlaps_r(inv, LEFT(tree), results));
    };
}
