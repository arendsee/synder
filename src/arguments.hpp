#ifndef __ARGUMENTS_H__
#define __ARGUMENTS_H__

#include <getopt.h>

#include "global.hpp"
#include "version.hpp"

#define MAX_POS 5

class Arguments
{

private:
    void set_defaults();
    void check_file(FILE * fp, char * name);
    void set_offsets(char * offset);

public:
    FILE *synfile;
    FILE *intfile;
    FILE *hitfile;
    std::string db_filename;
    std::string cmd;
    std::vector<std::string> pos;
    int offsets[4];
    long k;
    char trans;
    bool validate;
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
