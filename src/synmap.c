#include <assert.h>
#include <errno.h>

#include "synmap.h"

Synmap *init_Synmap()
{
  Synmap *syn = (Synmap *) malloc(sizeof(Synmap));
  syn->genome = (Genome **) malloc(2 * sizeof(Genome *));
  return (syn);
}

void free_Synmap(Synmap * synmap)
{
  if (synmap != NULL) {
    free_Genome(SG(synmap, 0));
    free_Genome(SG(synmap, 1));
    free(synmap->genome);
    free(synmap);
  }
}

void print_Synmap(Synmap * synmap, bool forward)
{
  print_Genome(SG(synmap, 0), forward);
  print_Genome(SG(synmap, 1), forward);
}

void sort_all_contigs(Synmap * synmap)
{
  for (int genid = 0; genid < 2; genid++) {
    for (int conid = 0; conid < SG(synmap, 0)->size; conid++) {
      sort_blocks_by_start(SGC(synmap, genid, conid));
      sort_blocks_by_stop(SGC(synmap, genid, conid));
    }
  }
}
