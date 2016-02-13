/** /todo figure out which of these I need */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "search.h"
#include "iv.h"

uint count_interval_overlaps_r(Interval *, struct IntervalTree *, uint);
uint count_point_overlaps_r(uint, struct IntervalTree *, uint);
IV * get_interval_overlaps_r(Interval *, struct IntervalTree *, IV *);
IV * get_point_overlaps_r(uint, struct IntervalTree *, IV *);

uint count_interval_overlaps(Interval * inv, struct IntervalTree * tree){
    return count_interval_overlaps_r(inv, tree, 0);
}

uint count_point_overlaps(uint point, struct IntervalTree * tree){
    return count_point_overlaps_r(point, tree, 0);
}

IA * get_point_overlaps(uint point, struct IntervalTree * tree){
    IV * results = get_point_overlaps_r(point, tree, iv_init(8));
    IA * ia = init_ia();
    /** \todo In itree, when I convert from iv to ia, the unused memory from iv
     * is hidden. I could free it be calling realloc.*/
    ia->size = results->size;
    ia->v = results->data;
    return ia;
}

IA * get_interval_overlaps(Interval * inv, struct IntervalTree * tree){
    IV * results = get_interval_overlaps_r(inv, tree, iv_init(8));
    IA * ia = init_ia();
    ia->size = results->size;
    ia->v = results->data;
    return ia;
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

IV * get_point_overlaps_r(uint point, struct IntervalTree * tree, IV * results){
    if(point >= tree->center) {
        for(int i = tree->by_stop->size - 1; i >= 0 ; i--){
            if(point <= tree->by_stop->v[i].stop){
                iv_add(results, tree->by_stop->v[i]);
            } else {
                break;
            }
        }
        if(tree->r_child)
            return get_point_overlaps_r(point, tree->r_child, results);
    }
    else {
        for(int i = 0; i < tree->by_start->size; i++){
            if(point >= tree->by_start->v[i].start){
                iv_add(results, tree->by_start->v[i]);
            } else {
                break;
            }
        }
        if(tree->l_child)
            return get_point_overlaps_r(point, tree->l_child, results);
    }
    return results;
}

IV * get_interval_overlaps_r(Interval * inv, struct IntervalTree * tree, IV * results){
    if(tree == NULL)
        return results;

    Pos center_location = point_overlap(tree->center, *inv);

    if(center_location == lo){
        for(int i = tree->by_stop->size - 1; i >= 0 ; i--){
            if(inv->start <= tree->by_stop->v[i].stop){
                iv_add(results, tree->by_stop->v[i]);
            } else {
                break;
            }
        }
        return get_interval_overlaps_r(inv, tree->r_child, results);
    }
    else if(center_location == hi){
        for(int i = 0; i < tree->by_start->size; i++){
            if(inv->stop >= tree->by_start->v[i].start){
                iv_add(results, tree->by_start->v[i]);
            } else {
                break;
            }
        }
        return get_interval_overlaps_r(inv, tree->l_child, results);
    }
    else{
        /** \todo In itree, make explicit function for splicing an IV and IA */
        for(int i = 0; i < tree->by_start->size; i++){
            iv_add(results, tree->by_start->v[i]);
        }
        return get_interval_overlaps_r(inv, tree->r_child,
               get_interval_overlaps_r(inv, tree->l_child, results));
    }
}
