#include "string.h"
#include "stdlib.h"
#include "stdio.h"

#include "genome.h"

Genome *init_genome(char *name, size_t size)
{
  Genome *gen = (Genome *) malloc(sizeof(Genome));
  gen->name = strdup(name);
  gen->size = size;
  gen->contig = (Contig **) malloc(size * sizeof(Contig *));
  return (gen);
}

void free_genome(Genome * genome)
{
  if (genome) {
    for (int i = 0; i < genome->size; i++) {
      free_contig(genome->contig[i]);
    }
    free(genome->contig);
    free(genome->name);
    free(genome);
  }
}

void print_genome(Genome * genome)
{
  printf(">\t%s\t%lu\n", genome->name, genome->size);
  for (int i = 0; i < genome->size; i++) {
    print_contig(genome->contig[i]);
  }
}
