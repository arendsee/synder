#ifndef __SIMULATE_H__
#define __SIMULATE_H__

#include <stdbool.h>

typedef uint unsigned int;

struct {
    uint start;
    uint stop;
    bool strand;
    size_t nlinks;
    struct S_Interval ** links;
} S_Interval;

struct {
    size_t size;
    struct S_Interval ** intervals;
} S_Contig

struct {
    size_t size;
    struct S_Contig ** contigs;
} S_Genome

struct {
    struct * S_Genome a;
    struct * S_Genome b;
} S_Synmap;

struct {
    // Likelihood of each event (should sum to 1)
    double flip_rate;
    double copy_rate;
    double move_rate;
    double delete_rate;
    double repeat_rate;

    // Length of affected region
    double flip_width;
    double copy_width;
    double move_width;
    double delete_width;
    double repeat_width;
} SimulationParameters;

struct SimulationParameters default_parameters();

struct S_Synmap init_s_synmap(uint ncontigs, double geom_rate);
struct S_Genome init_s_genome();
struct S_Contig init_s_contig();
struct S_Interval init_s_interval();

void free_s_synmap(struct S_Synmap * synmap);
void free_s_genome(struct S_Genome * genome);
void free_s_contig(struct S_Contig * contig);
void free_s_interval(struct S_Interval * interval);

void print_s_synmap(struct S_Synmap * synmap);
void print_s_genome(struct S_Genome * genome);
void print_s_contig(struct S_Contig * contig);
void print_s_interval(struct S_Interval * interval);

void s_mutate(struct S_Synmap * synmap, struct SimulationParameters * par);

void s_flip(struct S_Synmap * synmap, width);

void s_copy(struct S_Synmap * synmap, width);

void s_move(struct S_Synmap * synmap, width);

void s_delete(struct S_Synmap * synmap, width);

void s_repeat(struct S_Synmap * synmap, width);

#endif
