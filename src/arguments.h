#ifndef __ARGUMENTS_H__
#define __ARGUMENTS_H__

#include <getopt.h>
#include <array>

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
    void check_file(FILE* fp, char* name);
    void set_offsets(char* offset);

public:

    Command cmd = C_UNSET;
    FILE *synfile = nullptr;
    FILE *intfile = nullptr;
    FILE *hitfile = nullptr;
    FILE *tclfile = nullptr;
    FILE *qclfile = nullptr;
    bool dump_blks = false;
    bool swap = false;
    bool debug = false;
    std::array<int, 4> offsets {{ 0 }};
    long k = 0;
    char trans = 'i';
    std::vector<std::string> pos;

    Arguments(int argc, char *argv[]);
    ~Arguments();

    void print();
    void print_help();
    void print_version();

};

#endif
