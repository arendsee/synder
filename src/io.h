#ifndef __SYNMAP_IO_H__
#define __SYNMAP_IO_H__

#include "global.h"
#include "synmap.h"

/** Build synteny tree from specially formatted file.
 *
 * @warning This function is VERY picky about input. It expects input to be
 * formatted exactly as util/prepare-data.sh produces. You must not feed this
 * function raw synteny files. I currently have no input checks.
 *
 * @param synfile specially formatted synteny file
 *
 * @return pointer to a complete Synmap object
 */
Synmap *load_Synmap(FILE * synfile, int swap);

#endif
