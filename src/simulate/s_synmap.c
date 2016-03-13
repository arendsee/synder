#include "synmap.h"

struct S_Synmap init_s_synmap(uint ncontigs, double geom_rate){ }

void free_s_synmap(struct S_Synmap * synmap){ }

void print_s_synmap(struct S_Synmap * synmap){ }

void s_mutate(struct S_Synmap * synmap, struct SimulationParameters * par){ }

void s_flip(struct S_Synmap * synmap, width){ }

void s_copy(struct S_Synmap * synmap, width){ }

void s_move(struct S_Synmap * synmap, width){ }

void s_delete(struct S_Synmap * synmap, width){ }

void s_repeat(struct S_Synmap * synmap, width){ }
