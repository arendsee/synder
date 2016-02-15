/** /todo figure out which of these I need */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "search.h"
#include "iv.h"

uint count_interval_overlaps_r(Interval *, struct IntervalTree *, uint);
uint count_point_overlaps_r(uint, struct IntervalTree *, uint);

void get_interval_overlaps_r(Interval *, struct IntervalTree *, IntervalResult *);
void get_point_overlaps_r(uint, struct IntervalTree *, IntervalResult *);

uint count_interval_overlaps(Interval * inv, struct IntervalTree * tree){
    return count_interval_overlaps_r(inv, tree, 0);
}

uint count_point_overlaps(uint point, struct IntervalTree * tree){
    return count_point_overlaps_r(point, tree, 0);
}

IntervalResult * init_IntervalResult(size_t size){
    IntervalResult * res = (IntervalResult *)malloc(sizeof(IntervalResult));
    res->iv = iv_init(size);
    res->inbetween = false;
    return(res);
}

void free_IntervalResult(IntervalResult * res){
    if(res){
        if(res->iv)
            iv_free(res->iv);
        free(res);
    }
}

IntervalResult * get_point_overlaps(uint point, struct IntervalTree * tree){
    IntervalResult * res = init_IntervalResult(IV_INITIAL_SIZE);
    get_point_overlaps_r(point, tree, res);
    return res;
}

IntervalResult * get_interval_overlaps(Interval * inv, struct IntervalTree * tree){
    IntervalResult * res = init_IntervalResult(IV_INITIAL_SIZE);
    get_interval_overlaps_r(inv, tree, res);
    return res;
}



uint count_point_overlaps_r(uint point, struct IntervalTree * tree, uint count){
    if(point >= tree->center) {
        for(int i = tree->by_stop->size - 1; i >= 0 ; i--){
            if(point <= tree->by_stop->v[i].stop){
                count++;
            } else {
                break;
            }
        }
        if(tree->r_child)
            return count_point_overlaps_r(point, tree->r_child, count);
    }
    else {
        for(int i = 0; i < tree->by_start->size; i++){
            if(point >= tree->by_start->v[i].start){
                count++;
            } else {
                break;
            }
        }
        if(tree->l_child)
            return count_point_overlaps_r(point, tree->l_child, count);
    }
    return count;
}

uint count_interval_overlaps_r(Interval * inv, struct IntervalTree * tree, uint count){
    if(tree == NULL)
        return count;

    Pos center_location = point_overlap(tree->center, *inv);

    if(center_location == lo){
        for(int i = tree->by_stop->size - 1; i >= 0 ; i--){
            if(inv->start <= tree->by_stop->v[i].stop){
                count++;
            } else {
                break;
            }
        }
        return count_interval_overlaps_r(inv, tree->r_child, count);
    }
    else if(center_location == hi){
        for(int i = 0; i < tree->by_start->size; i++){
            if(inv->stop >= tree->by_start->v[i].start){
                count++;
            } else {
                break;
            }
        }
        return count_interval_overlaps_r(inv, tree->l_child, count);
    }
    else{
        count += tree->by_start->size;
        return count_interval_overlaps_r(inv, tree->r_child,
               count_interval_overlaps_r(inv, tree->l_child, count));
    }
}

void get_point_overlaps_r(uint point, struct IntervalTree * tree, IntervalResult * results){
    if(point >= tree->center) {
        for(int i = tree->by_stop->size - 1; i >= 0 ; i--){
            if(point <= tree->by_stop->v[i].stop){
                iv_add(results->iv, tree->by_stop->v[i]);
            } else {
                break;
            }
        }
        if(tree->r_child){
            get_point_overlaps_r(point, tree->r_child, results);
        } 
        else if(!results->iv->size) {
            results->inbetween = true;
            iv_add(results->iv, LAST_STOP(tree));
        }
    }
    else {
        for(int i = 0; i < tree->by_start->size; i++){
            if(point >= tree->by_start->v[i].start){
                iv_add(results->iv, tree->by_start->v[i]);
            } else {
                break;
            }
        }
        if(tree->l_child)
            get_point_overlaps_r(point, tree->l_child, results);
        else if(!results->iv->size) {
            results->inbetween = true;
            iv_add(results->iv, FIRST_START(tree));
        }
    }
}

void get_interval_overlaps_r(Interval * inv, struct IntervalTree * tree, IntervalResult * results){
    Pos center_location = point_overlap(tree->center, *inv);

    if(center_location == lo){
        for(int i = tree->by_stop->size - 1; i >= 0 ; i--){
            if(inv->start <= tree->by_stop->v[i].stop){
                iv_add(results->iv, tree->by_stop->v[i]);
            } else {
                break;
            }
        }
        if(tree->r_child){
            return get_interval_overlaps_r(inv, tree->r_child, results);
        } else if(!results->iv->size) {
            results->inbetween = true;
            iv_add(results->iv, LAST_STOP(tree));
        }
    }
    else if(center_location == hi){
        for(int i = 0; i < tree->by_start->size; i++){
            if(inv->stop >= tree->by_start->v[i].start){
                iv_add(results->iv, tree->by_start->v[i]);
            } else {
                break;
            }
        }
        if(tree->l_child){
            get_interval_overlaps_r(inv, tree->l_child, results);
        } else if(!results->iv->size){
            results->inbetween = true;
            iv_add(results->iv, FIRST_START(tree));
        }
    }
    else{
        /** \todo In itree, make explicit function for splicing an IV and IA */
        for(int i = 0; i < tree->by_start->size; i++){
            iv_add(results->iv, tree->by_start->v[i]);
        }
        get_interval_overlaps_r(inv, tree->r_child, results);
        get_interval_overlaps_r(inv, tree->l_child, results);
    }
}
