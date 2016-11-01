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


//' print all blocks with contiguous set ids
//'
//' @param filename synteny map file name
// [[Rcpp::export]]
Rcpp::DataFrame c_dump (std::string filename)
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

    return syn.as_data_frame();
}

//' predict search intervals
//'
//' @param synfilename synteny map file name
//' @param gfffilename GFF file name
// [[Rcpp::export]]
Rcpp::DataFrame c_search(std::string synfilename, std::string gfffilename)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* gfffile = fopen(gfffilename.c_str(), "r");

    Synmap syn(
        synfile, // FILE* synfile
        nullptr, // FILE* tclfile
        nullptr, // FILE* qclfile
        false,   // bool swap
        0,       // int k
        0.001,   // double r
        'i'      // char trans
    );

    return syn.search(gfffile);
}


//' remove links that disagree with the synteny map
//'
//' @param synfilename synteny map file name
//' @param intfilename int file name
// [[Rcpp::export]]
Rcpp::CharacterVector c_filter(std::string synfilename, std::string intfilename)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* intfile = fopen(intfilename.c_str(), "r");

    Synmap syn(
        synfile, // FILE* synfile
        nullptr, // FILE* tclfile
        nullptr, // FILE* qclfile
        false,   // bool swap
        0,       // int k
        0.001,   // double r
        'i'      // char trans
    );

    return syn.filter(intfile);
}

//' trace intervals across genomes
//'
//' @param synfilename synteny map file name
//' @param gfffilename gff file name
// [[Rcpp::export]]
Rcpp::DataFrame c_map(std::string synfilename, std::string gfffilename)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* gfffile = fopen(gfffilename.c_str(), "r");

    Synmap syn(
        synfile, // FILE* synfile
        nullptr, // FILE* tclfile
        nullptr, // FILE* qclfile
        false,   // bool swap
        0,       // int k
        0.001,   // double r
        'i'      // char trans
    );

    return syn.map(gfffile);
}

//' count overlaps
//'
//' @param synfilename synteny map file name
//' @param gfffilename gff file name
// [[Rcpp::export]]
Rcpp::DataFrame c_count(std::string synfilename, std::string gfffilename)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* gfffile = fopen(gfffilename.c_str(), "r");

    Synmap syn(
        synfile, // FILE* synfile
        nullptr, // FILE* tclfile
        nullptr, // FILE* qclfile
        false,   // bool swap
        0,       // int k
        0.001,   // double r
        'i'      // char trans
    );

    return syn.count(gfffile);
}
