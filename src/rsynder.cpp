#include <string>
#include <iostream>

#include "global.h"
#include "synmap.h"

int Offsets::syn_start;
int Offsets::syn_stop;
int Offsets::in_start;
int Offsets::in_stop;
int Offsets::out_start;
int Offsets::out_stop;


//' Dump the synteny map
//'
//' @param filename synteny map file name
//' @export
// [[Rcpp::export]]
void dump (std::string filename)
{
    FILE* synfile = fopen(filename.c_str(), "r");

    Synmap syn(
        synfile, // FILE* synfile
        nullptr, // FILE* tclfile
        nullptr, // FILE* qclfile
        false,   // bool swap
        0,       // int k
        0.001,   // double r
        'i'      // char trans
    );
    syn.dump_blocks();
}



    // "  -i, --input \tGFF file of input intervals (default: STDIN)"
    // "  -i, --input \tA TAB-delimited map file (default: STDIN)\n"
    // "              \tThe first 6 fields must be:\n"
    // "              \t  1. qseqid - query contig id (e.g. Chr1)\n"
    // "              \t  2. qstart - query interval start\n"
    // "              \t  3. qstop  - query interval stop\n"
    // "              \t  4. sseqid - target contig id\n"
    // "              \t  5. sstart - target interval start\n"
    // "              \t  6. sstop  - target interval stop\n"
    // "              \t"
    // "Any number of additional columns may follow, they will be printed without alteration"
    // "  -s, --synmap \tThe primary synteny map (e.g. output of Satsuma)\n"
    // "               \tTAB-delimited, no header, with fields\n"
    // "               \t  1. qseqid - query contig id (e.g. Chr1)\n"
    // "               \t  2. qstart - query interval start\n"
    // "               \t  3. qstop  - query interval stop\n"
    // "               \t  4. sseqid - target contig id\n"
    // "               \t  5. sstart - target interval start\n"
    // "               \t  6. sstop  - target interval stop\n"
    // "               \t  7. score  - score of the syntenic match*\n"
    // "               \t  8. strand - relative orientation\n"
    // "               \t  * score can be any numeric value, it will be\n"
    // "               \t    transformed as specified by the -x option"
    // "  -t, --tcl \tTarget chromosome lengths file\n"
    // "            \tTAB-delimited with columns: <name>, <length>"
    // "  -q, --qcl \tQuery chromosome lengths file\n"
    // "            \tTAB-delimited with columns: <name>, <length>"
    // "  -r, --reverse \tSwitch target and query in synteny map"
    // "  -b, --base-offsets \tStart and stop offsets for synteny map, input and output (e.g. 011100)"
    // "  -k, --interruption-fuzz \tNumber of interrupting intervals allowed before breaking contiguous set (default=0)"
    // "  -j, --decay-rate \tScore decay rate (default=0.001)"
    // "  -x, --transform \tTransform score (Synder requires additive scores):\n"
    // "                 \t -'i' := S            (default, no transformation)\n"
    // "                 \t -'d' := L * S        (transform from score densities)\n"
    // "                 \t -'p' := L * S / 100  (transform from percent identity)\n"
    // "                 \t -'l' := -log(S)      (transform from e-values or p-values)\n"
    // "                 \t   Where S is input score and L interval length"
    //
    //         "synder dump - print all blocks with contiguous set ids\n"
    //         "USAGE\n"
    //         "  synder dump -s SYNTENY_MAP\n"
    //         "DESCRIPTION\n"
    //         "  Dump processed blocks to stdout and exit. The output contains\n"
    //         "  all the fields in the synteny map, with preserved order, and\n"
    //         "  adds a column of contiguous set ids. Also, doubly overlapping\n"
    //         "  bocks are also merged.\n"
    //         "ARGUMENTS"
    //
    //         "synder filter - remove links that disagree with the synteny map\n"
    //         "USAGE\n"
    //         "  synder filter -i INPUT.gff -s SYNMAP.tab [OPTIONS]\n"
    //         "DESCRIPTION\n"
    //         "  Remove input links that disagree with the synteny map Given a map\n"
    //         "  between the target and query (e.g. BLAST output) and a synteny map,\n"
    //         "  remove all entries in the map that do not overlap search intervals\n"
    //         "  in the query.\n"
    //         "ARGUMENTS"
    //
    //         "synder map - trace intervals across genomes\n"
    //         "USAGE\n"
    //         "  synder map -i INPUT.gff -s SYNMAP.tab [OPTIONS]\n"
    //         "DESCRIPTION\n"
    //         "  Given a set of query intervals, find all synmap queries that\n"
    //         "  overlap and map to homologous target intervals.\n"
    //         "ARGUMENTS"
    //
    //         "synder count - count overlaps\n"
    //         "USAGE\n"
    //         "  synder count -i INPUT.gff -s SYNMAP.tab [OPTIONS]\n"
    //         "DESCRIPTION\n"
    //         "  Given a set of query intervals, count all synmap queries that\n"
    //         "  overlap and map to homologous target intervals.\n"
    //         "ARGUMENTS"
    //
    //         "synder search - predict search intervals\n"
    //         "USAGE\n"
    //         "  synder search -i INPUT.gff -s SYNMAP.tab [OPTIONS]\n"
    //         "DESCRIPTION\n"
    //         "  Given an input GFF and a synteny, find search intervals. See README for\n"
    //         "  algorithmic details and a discussion of scoring and argument -k.\n"
    //         "ARGUMENTS"
    //         "EXAMPLE\n"
    //         "  # Using files in the sample-inputs/ folder\n"
    //         "  synder search -i a.gff -s a-b.tab -b 0011 -x p\n"
    //         "  # Search the opposite direction\n"
    //         "  synder search -r -i b.gff -s a-b.tab\n"


// void set_arguments(Arguments& arg, std::vector<option::Option> options)
// {
//     if(options[INPUT_GFF].desc != 0) {
//         arg.intfile = check_open(options[INPUT_GFF].arg, "GFF file");
//     }
//
//     if(options[INPUT_MAP].desc != 0) {
//         arg.intfile = check_open(options[INPUT_MAP].arg, "input map");
//     }
//
//     if(options[SYNMAP].desc != 0) {
//         arg.synfile = check_open(options[SYNMAP].arg, "synteny map");
//     }
//
//     if(options[TCLFILE].desc != 0) {
//         arg.tclfile = check_open(options[TCLFILE].arg, "target contig length file");
//     }
//
//     if(options[QCLFILE].desc != 0) {
//         arg.qclfile = check_open(options[QCLFILE].arg, "query contig length file");
//     }
//
//     if(options[REVERSE].desc != 0) {
//         arg.swap = true;
//     }
//
//     if(options[BASE_OFFSETS].desc != 0) {
//         std::string a(options[BASE_OFFSETS].arg);
//         if(a.length() == 6) {
//             Offsets::syn_start  = a[0] == '0' ? 0 : 1;
//             Offsets::syn_stop   = a[1] == '0' ? 0 : 1;
//             Offsets::in_start   = a[2] == '0' ? 0 : 1;
//             Offsets::in_stop    = a[3] == '0' ? 0 : 1;
//             Offsets::out_start  = a[4] == '0' ? 0 : 1;
//             Offsets::out_stop   = a[5] == '0' ? 0 : 1;
//         } else {
//             fprintf(stderr, "Offsets must be a string of length 6\n");
//             exit(EXIT_FAILURE);
//         }
//     }
//
//     if(options[K].desc != 0) {
//         long k = strtol(options[K].arg, nullptr, 0);
//         if ((errno == ERANGE && (k == LONG_MAX || k == LONG_MIN)) || (errno != 0 && k == 0)) {
//             perror("In argument -k NUM, NUM must be an integer");
//         }
//         arg.k = k;
//     }
//
//     if(options[TRANSFORM].desc != 0) {
//         char trans = options[TRANSFORM].arg[0];
//         if (! (trans == 'i' || trans == 'd' || trans == 'p' || trans == 'l')) {
//             fprintf(stderr, "-x only takes arguments 'i', 'd', 'p' and 'l'\n");
//             exit(EXIT_FAILURE);
//         }
//         arg.trans = trans;
//     }
//
//     // If no file given by -i, use STDIN
//     if (arg.intfile == nullptr) {
//         arg.intfile = stdin;
//     }
// }
//
// FILE* check_open(const char* name, const char* adjective)
// {
//     FILE* fp = fopen(name, "r");
//     if (fp == nullptr) {
//         fprintf(stderr, "ERROR: Failed to open %s '%s'\n", adjective, name);
//         exit(EXIT_FAILURE);
//     }
//     return fp;
// }
