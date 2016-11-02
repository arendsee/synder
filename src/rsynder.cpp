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

// [[Rcpp::export]]
Rcpp::IntegerVector c_logical_strand (Rcpp::CharacterVector cv) {
    size_t N = cv.size();
    Rcpp::IntegerVector iv(N);
    for(int i = 0; i < N; i++){
        if(cv[i] == "+"){
            iv[i] = 1;
        } else if (cv[i] == "-"){
            iv[i] = 0;
        } else {
            iv[i] = NA_INTEGER;
        }
    }
    // Recast as a logical vector in the R-side function
    return iv;
}

//' print all blocks with contiguous set ids
//'
//' @param filename synteny map file name
//' @param swap     reverse direction of synteny map (e.g. swap query and target) 
//' @param trans    STUB
// [[Rcpp::export]]
Rcpp::DataFrame c_dump (
    std::string filename,
    bool swap,
    char trans
)
{
    FILE* synfh = fopen(filename.c_str(), "r");

    Synmap syn(synfh, nullptr, nullptr, swap, 0, 0, trans);

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
    char trans
)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* gfffile = fopen(gfffilename.c_str(), "r");
    FILE* tclfile = fopen(tclfilename.c_str(), "r");
    FILE* qclfile = fopen(qclfilename.c_str(), "r");

    Synmap syn(synfile, tclfile, qclfile, swap, k, r, trans);

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
    char trans
)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* intfile = fopen(intfilename.c_str(), "r");

    Synmap syn(synfile, nullptr, nullptr, swap, k, r, trans);

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
    bool swap
)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* gfffile = fopen(gfffilename.c_str(), "r");

    Synmap syn(synfile, nullptr, nullptr, swap, 0, 0, 'i');

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
    bool swap
)
{
    FILE* synfile = fopen(synfilename.c_str(), "r");
    FILE* gfffile = fopen(gfffilename.c_str(), "r");

    Synmap syn(synfile, nullptr, nullptr, swap, 0, 0, 'i');

    return syn.count(gfffile);
}
