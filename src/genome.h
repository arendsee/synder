#ifndef __GENOME_H__
#define __GENOME_H__

#include <stdbool.h>

#include "contig.h"

/** A named set of Contig objects*/
typedef struct {
  char *name;
  size_t size;
  Contig **contig;
} Genome;

/** Allocate memory for Genome *name* of size *size*.
 *
 * @param name genome name (e.g. "Arabidopsis_thaliana")
 * @param size number of child Contig structs (e.g. chromosomes or scaffolds) 
 *
 * @return pointer to new Genome struct 
 * */
Genome *init_genome(char *, size_t);

/**
 * Recursively free all memory
 *
 * For each Contig in the contig field, calls free_contig.
 *
 * @param pointer to a Genome struct
 */
void free_genome(Genome *);

/** Recursively print a genome. */
void print_genome(Genome *, bool forward);

#endif
