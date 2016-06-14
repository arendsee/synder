#define _GNU_SOURCE
#include <getopt.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


#include "ui.h"

#define MAX_POS 3

/** contains all legal arguments from all subcommands */
Arguments create_Arguments() {
    Arguments args = {
        .synfile = NULL,
        .intfile = NULL,
        .hitfile = NULL,
        .db_filename = NULL,
        .cmd = NULL,
        .pos = (char **)malloc(MAX_POS * sizeof(char *)),
        .test = false
    };
    memset(args.pos, 0, MAX_POS * sizeof(char *));
    return args;
}

void close_Arguments(Arguments arg){
    if(arg.synfile)
        fclose(arg.synfile);
    if(arg.intfile)
        fclose(arg.intfile);
    if(arg.hitfile)
        fclose(arg.hitfile);
    if(arg.db_filename)
        free(arg.db_filename);
    if(arg.pos){
        for(int i = 0; i < 3; i++){
            if(arg.pos[i])
                free(arg.pos[i]);
        }
        free(arg.pos);
    }
    if(arg.cmd)
        free(arg.cmd);
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
    printf(
    "USAGE\n"
    "  synder -d SYNTENY_FILE QUERY TARGET DB_DIR\n"
    "  synder [OPTIONS] -s SYNTENY_DB -c COMMAND\n"
    "COMMANDS\n"
    "  map    - print target intervals overlapping each query interval\n"
    "  count  - like map but prints only the number that overlap\n"
    "  search - find target search spaces for each query interval\n"
    "  filter - print query-to-target links consistent with the synteny map\n"
    "EXAMPLES\n"
    "  synder -d at-al.tab at al db\n"
    "  synder -i at.gff   -s db/at_al.txt -c map\n"
    "  synder -i at.gff   -s db/at_al.txt -c count\n"
    "  synder -i at.syn   -s db/at_al.txt -c search\n"
    "  synder -f hits.syn -s db/at_al.txt -c filter\n"
    "  synder test\n"
    );
    exit(EXIT_SUCCESS);
}

Arguments parse_command(int argc, char * argv[]){
    if(argc == 1)
        print_help();
    int opt;
    FILE * temp;
    Arguments args = create_Arguments();
    if(argc == 2 && strcmp(argv[1], "test") == 0){
        args.test = true; 
        return args;
    }
    while((opt = getopt(argc, argv, "hd:s:i:c:f:")) != -1){
        switch(opt) {
            case 'h':
                print_help();
                break;
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
            case 'f':
                args.hitfile = fopen(optarg, "r");
                check_file(args.hitfile, optarg);
                break;
            case 'c':
                args.cmd = strdup(optarg);
                break;
            case '?':
                exit(EXIT_FAILURE);
        }
    }

    for(int i = 0; optind < argc; optind++, i++){		
        if(i >= MAX_POS){		
            printf("There can be only %d positionals\n", MAX_POS);		
            exit(EXIT_FAILURE);		
        }		
        args.pos[i] = strdup(argv[optind]);		
    }
    return args;
}
