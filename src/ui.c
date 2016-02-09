#define _GNU_SOURCE
#include <getopt.h>

#include "ui.h"

Arguments create_Arguments() {
    Arguments args = {
        .chr   = 0,
        .start = 0,
        .stop  = 0,
        .synfile = NULL
    };
    return args;
}

void close_Arguments(Arguments arg){
    if(arg.synfile)
        fclose(arg.synfile);
}

void print_args(Arguments args){
    printf("chr-%lu (%lu, %lu)\n", args.chr, args.start, args.stop);
}

Arguments parse_command(int argc, char * argv[]){
    int opt;
    Arguments args = create_Arguments();
    while((opt = getopt(argc, argv, "a:b:c:f:")) != -1){
        switch(opt) {
            case 'a':
                args.chr = atoi(optarg);
                break;
            case 'b': 
                args.start = atoi(optarg);
                break;
            case 'c': 
                args.stop = atoi(optarg);
                break;
            case 'f':
                args.synfile = fopen(optarg, "r");
                break;
            case '?':
                exit(EXIT_FAILURE);
        }
    }

    if(args.synfile == NULL){
        fprintf(stderr, "ERROR: Failed to open synteny file '%s'\n", optarg);
        exit(EXIT_FAILURE);
    }

    for(; optind < argc; optind++){
        printf("Positional: %s\n", argv[optind]);
    }
    return args;
}
