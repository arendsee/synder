#include "string.h"
#include "stdlib.h"
#include "stdio.h"

#include "genome.h"

Genome *init_Genome(char *name, size_t size)
{
  Genome *gen = (Genome *) malloc(sizeof(Genome));
  gen->name = strdup(name);
  gen->size = size;
  gen->contig = (Contig **) malloc(size * sizeof(Contig *));
  return (gen);
}

void free_Genome(Genome * genome)
{
  if (genome != NULL) {
    for (int i = 0; i < genome->size; i++) {
      free_Contig(genome->contig[i]);
    }
    free(genome->contig);
    free(genome->name);
    free(genome);
  }
}

void print_Genome(Genome * genome, bool forward)
{
  printf(">\t%s\t%lu\n", genome->name, genome->size);
  for (int i = 0; i < genome->size; i++) {
    print_Contig(genome->contig[i], forward);
  }
}
