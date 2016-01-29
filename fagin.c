#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define _GNU_SOURCE
#include <getopt.h>

typedef struct Arguments {
    unsigned int a:1;
    unsigned int b:1;
    unsigned int c:1;
    char * file;
} Arguments;

Arguments create_Arguments() {
    Arguments args = {
        .a = false,
        .b = false,
        .c = false
    };
    return args;
}

int main(int argc, char * argv[]){
    int opt;
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
                strcpy(args.file, optarg);
                break;
            case '?':
                exit(EXIT_FAILURE);
        }
    }

    for(; optind < argc; optind++){
        printf("Positional: %s\n", argv[optind]);
    }

    return(EXIT_SUCCESS);
}
