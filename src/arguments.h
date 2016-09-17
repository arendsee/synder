#ifndef __ARGUMENTS_H__
#define __ARGUMENTS_H__

#include <getopt.h>

#include "global.h"
#include "version.h"

#define MAX_POS 5

typedef enum {
    C_UNSET,
    C_FILTER,
    C_DEBUG,
    C_DUMP,
    C_MAP,
    C_COUNT,
    C_SEARCH
} Command;

class Arguments
{

private:
    void set_defaults();
    void check_file(FILE* fp, char* name);
    void set_offsets(char* offset);

public:
    FILE *synfile;
    FILE *intfile;
    FILE *hitfile;
    FILE *tclfile;
    FILE *qclfile;
    Command cmd;
    std::vector<std::string> pos;
    int offsets[4];
    long k;
    char trans;
    bool dump_blks;
    bool swap;
    bool debug;

    Arguments(int argc, char *argv[]);
    ~Arguments();
    void print();
    void print_help();
    void print_version();

};

#endif
