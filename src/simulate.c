struct SimulationParameters default_parameters(){ }

struct S_Synmap init_s_synmap(uint ncontigs, double geom_rate){ }

struct S_Genome init_s_genome(){ }

struct S_Contig init_s_contig(){ }

struct S_Interval init_s_interval(){ }

void free_s_synmap(struct S_Synmap * synmap){ }

void free_s_genome(struct S_Genome * genome){ }

void free_s_contig(struct S_Contig * contig){ }

void free_s_interval(struct S_Interval * interval){ }

void print_s_synmap(struct S_Synmap * synmap){ }

void print_s_genome(struct S_Genome * genome){ }

void print_s_contig(struct S_Contig * contig){ }

void print_s_interval(struct S_Interval * interval){ }

void s_mutate(struct S_Synmap * synmap, struct SimulationParameters * par){ }

void s_flip(struct S_Synmap * synmap, width){ }

void s_copy(struct S_Synmap * synmap, width){ }

void s_move(struct S_Synmap * synmap, width){ }

void s_delete(struct S_Synmap * synmap, width){ }

void s_repeat(struct S_Synmap * synmap, width){ }
