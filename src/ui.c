#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define _GNU_SOURCE
#include <getopt.h>

#include "ui.h"

Arguments create_Arguments() {
    Arguments args = {
        .a = false,
        .b = false,
        .c = false
    };
    return args;
}

void free_Arguments(Arguments arg){
    if(arg.file)
        free(arg.file);
}

void print_args(Arguments args){
    printf("a=%d b=%d c=%d file=%s\n",
            args.a ? 1 : 0,
            args.b ? 1 : 0,
            args.c ? 1 : 0,
            args.file);
}

Arguments parse_command(int argc, char * argv[]){
    int opt;
    size_t size;
    Arguments args = create_Arguments();
    while((opt = getopt(argc, argv, "abcf:")) != -1){
        switch(opt) {
            case 'a':
                args.a = true;
                break;
            case 'b': 
                args.b = true;
                break;
            case 'c': 
                args.c = true;
                break;
            case 'f':
                size = sizeof(char) * strlen(optarg);
                args.file = (char *)malloc(size + sizeof(char));
                strncpy(args.file, optarg, size);
                break;
            case '?':
                exit(EXIT_FAILURE);
        }
    }
    for(; optind < argc; optind++){
        printf("Positional: %s\n", argv[optind]);
    }
    return args;
}
