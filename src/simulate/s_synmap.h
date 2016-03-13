#ifndef __SYNMAP_H__
#define __SYNMAP_H__

struct {
    struct * S_Genome a;
    struct * S_Genome b;
} S_Synmap;

void free_s_synmap(struct S_Synmap * synmap);
void print_s_synmap(struct S_Synmap * synmap);

struct S_Synmap init_s_synmap(uint ncontigs, double geom_rate);

void s_mutate(struct S_Synmap * synmap, struct SimulationParameters * par);

void s_flip(struct S_Synmap * synmap, width);

void s_copy(struct S_Synmap * synmap, width);

void s_move(struct S_Synmap * synmap, width);

void s_delete(struct S_Synmap * synmap, width);

void s_repeat(struct S_Synmap * synmap, width);

#endif
