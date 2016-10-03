#ifndef __FEATURE_H__
#define __FEATURE_H__

#include "interval.hpp"
#include "global.h"

#include <string>
#include <iostream>
#include <climits>


class Feature : public Interval<Feature>
{
public:

    std::string parent_name = ".";
    std::string name        = ".";
    long parent_length      = DEFAULT_CONTIG_LENGTH;
    char strand             = '.';

    Feature() { }

    Feature(
        const char* t_parent_name,
        long        t_start,
        long        t_stop
    )
        :
        Interval(t_start, t_stop),
        parent_name(t_parent_name)
    { }

    Feature(
        const char* t_parent_name,
        long        t_start,
        long        t_stop,
        const char* t_name,
        long        t_parent_length,
        char        t_strand='+'
    )
        :
        Interval(t_start, t_stop),
        parent_name(t_parent_name),
        name(t_name),
        parent_length(t_parent_length),
        strand(t_strand)
    { }

    ~Feature() { };

    bool feature_overlap(Feature* other)
    {
        return overlap(other) && parent_name == other->parent_name;
    }

};

#endif
