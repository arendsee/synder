#define _GNU_SOURCE
#include <getopt.h>
#include <string.h>

#include "ui.h"

/** contains all legal arguments from all subcommands */
Arguments create_Arguments() {
    Arguments args = {
        .synfile = NULL,
        .intfile = NULL,
        .db_filename  = NULL
    };
    return args;
}

void close_Arguments(Arguments arg){
    if(arg.synfile)
        fclose(arg.synfile);
    if(arg.intfile)
        fclose(arg.intfile);
    if(arg.db_filename)
        free(arg.db_filename);
}

void print_args(Arguments args){
    printf("stub\n");
}

void check_file(FILE * fp, char * name){
    if(fp == NULL){
        fprintf(stderr, "ERROR: Failed to open '%s'\n", name);
        exit(EXIT_FAILURE);
    }
}

void print_help(){
    printf("USAGE: -d SYNTENY_FILE\n");
    printf("USAGE: -i GFF_FILE -s SYNTENY_DB\n");
}

Arguments parse_command(int argc, char * argv[]){
    int opt;
    FILE * temp;
    Arguments args = create_Arguments();
    while((opt = getopt(argc, argv, "hd:s:i:")) != -1){
        switch(opt) {
            case 'h':
                print_help();
                exit(EXIT_SUCCESS);
            case 'd':
                temp = fopen(optarg, "r");
                check_file(temp, optarg);
                fclose(temp);
                args.db_filename = strdup(optarg);
                break;
            case 's':
                args.synfile = fopen(optarg, "r");
                check_file(args.synfile, optarg);
                break;
            case 'i':
                args.intfile = fopen(optarg, "r");
                check_file(args.intfile, optarg);
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
