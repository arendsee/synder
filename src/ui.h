#ifndef __UI_H__
#define __UI_U__

#include <stdio.h>

typedef struct Arguments {
    unsigned int a:1;
    unsigned int b:1;
    unsigned int c:1;
    FILE * synfile;
} Arguments;

void close_Arguments(Arguments arg);

Arguments create_Arguments();

void print_args(Arguments args);

Arguments parse_command(int argc, char * argv[]);

#endif
