#include <string>
#include <iostream>

#include "global.h"
#include "synmap.h"

//' print all blocks with contiguous set ids
//'
//' @param filename synteny map file name
//' @param swap     reverse direction of synteny map (e.g. swap query and target) 
//' @param trans    STUB
// [[Rcpp::export]]
Rcpp::DataFrame c_dump (
    std::string filename,
    bool swap,
    char trans,
    std::vector<int> offsets
)
{
    FILE* synfh = fopen(filename.c_str(), "r");

    Synmap syn(synfh, nullptr, nullptr, swap, 0, 0, trans, offsets);

    return syn.as_data_frame();
}

//' predict search intervals
//'
//' @param synfilename synteny map file name
//' @param gfffilename GFF file name
//' @param swap     reverse direction of synteny map (e.g. swap query and target) 
//' @param k        STUB
//' @param r        STUB
//' @param trans    stuB
// [[Rcpp::export]]
Rcpp::DataFrame c_search(
    std::string synfilename,
    std::string gfffilename,
    std::string tclfilename,
    std::string qclfilename,
    bool swap,
    int k,
    double r,
    char trans,
    std::vector<int> offsets
)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* gfffile = fopen(gfffilename.c_str(), "r");
    FILE* tclfile = fopen(tclfilename.c_str(), "r");
    FILE* qclfile = fopen(qclfilename.c_str(), "r");

    Synmap syn(synfile, tclfile, qclfile, swap, k, r, trans, offsets);

    return syn.search(gfffile);
}


//' remove links that disagree with the synteny map
//'
//' @param synfilename synteny map file name
//' @param intfilename int file name
//' @param swap     reverse direction of synteny map (e.g. swap query and target) 
//' @param k        STUB
//' @param r        STUB
//' @param trans    stuB
// [[Rcpp::export]]
Rcpp::CharacterVector c_filter(
    std::string synfilename,
    std::string intfilename,
    bool swap,
    int k,
    double r,
    char trans,
    std::vector<int> offsets
)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* intfile = fopen(intfilename.c_str(), "r");

    Synmap syn(synfile, nullptr, nullptr, swap, k, r, trans, offsets);

    return syn.filter(intfile);
}

//' trace intervals across genomes
//'
//' @param synfilename synteny map file name
//' @param gfffilename gff file name
//' @param swap        reverse direction of synteny map (e.g. swap query and target) 
// [[Rcpp::export]]
Rcpp::DataFrame c_map(
    std::string synfilename,
    std::string gfffilename,
    bool swap,
    std::vector<int> offsets
)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* gfffile = fopen(gfffilename.c_str(), "r");

    Synmap syn(synfile, nullptr, nullptr, swap, 0, 0, 'i', offsets);

    return syn.map(gfffile);
}

//' count overlaps
//'
//' @param synfilename synteny map file name
//' @param gfffilename gff file name
//' @param swap        reverse direction of synteny map (e.g. swap query and target) 
// [[Rcpp::export]]
Rcpp::DataFrame c_count(
    std::string synfilename,
    std::string gfffilename,
    bool swap,
    std::vector<int> offsets
)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* gfffile = fopen(gfffilename.c_str(), "r");

    Synmap syn(synfile, nullptr, nullptr, swap, 0, 0, 'i', offsets);

    return syn.count(gfffile);
}
