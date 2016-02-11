#include <stdio.h>

#include "synmap.h"

typedef enum {COUNT, MAP} Command;

void analysis_count(Synmap * syn, FILE * intfile);

void analysis_map(Synmap * syn, FILE * intfile);
