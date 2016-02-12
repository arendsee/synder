#ifndef __SYNMAP_H__
#define __SYNMAP_H__

#include <stdio.h>

#include "genome/genome.h"

#ifndef uint
#define uint unsigned int
#endif

#define SGCB(syn, gid, cid, bid) syn->genome[(gid)]->contig[(cid)]->block[(bid)]
#define SGC(syn, gid, cid)       syn->genome[(gid)]->contig[(cid)]
#define SG(syn, gid)             syn->genome[(gid)]

#define QT_SGCB(syn, blk) syn->genome[1]->contig[blk->oseqid]->block[blk->oblkid]
#define QT_SGC(syn, blk)  syn->genome[1]->contig[blk->oseqid]

/** A pair of syntenically linked Genome objects  */
typedef struct {
    Genome ** genome;
} Synmap;

/**
 * Allocate memory for a new Synmap.
 *
 * This struct holds the pair of genomes that wil be compared.
 *
 * @return pointer to the new Synmap
 */
Synmap * init_synmap();

/** Recursively free all memory allocated to the synteny map.
 * 
 * Calls free_genome on both its Genome children.
 *
 * @param pointer to Synmap struct
 *
 * */
void free_synmap(Synmap *);

/** Recursively print a synteny map. */
void print_synmap(Synmap *);

#endif
