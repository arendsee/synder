#include <string>
#include <iostream>

#include "global.h"
#include "synmap.h"

//' print all blocks with contiguous set ids
//'
//' @param syn      synteny map file name
//' @param swap     reverse direction of synteny map (e.g. swap query and target) 
//' @param trans    score transform methods, single character
//' @param offsets  6-element integer vector of [01] offsets
// [[Rcpp::export]]
Rcpp::DataFrame c_dump (
    std::string syn,
    bool swap,
    char trans,
    std::vector<int> offsets
)
{
    FILE* synfh = fopen(syn.c_str(), "r");

    Synmap synmap(synfh, nullptr, nullptr, swap, 0, 0, trans, offsets);

    return synmap.as_data_frame();
}

//' predict search intervals
//'
//' @param syn     synteny map file name
//' @param gff     GFF file name
//' @param tcl     target contig lengths file name
//' @param qcl     target contig lengths file name
//' @param swap    reverse direction of synteny map (e.g. swap query and target) 
//' @param k       match fuziness, integer
//' @param r       score decay rate, 0 means no context, high means more context
//' @param trans   score transform methods, single character
//' @param offsets 6-element integer vector of [01] offsets
// [[Rcpp::export]]
Rcpp::DataFrame c_search(
    std::string syn,
    std::string gff,
    std::string tcl,
    std::string qcl,
    bool swap,
    int k,
    double r,
    char trans,
    std::vector<int> offsets
)
{
    FILE* synfh = fopen(syn.c_str(), "r");
    FILE* gfffh = fopen(gff.c_str(), "r");
    FILE* tclfh = fopen(tcl.c_str(), "r");
    FILE* qclfh = fopen(qcl.c_str(), "r");

    Synmap synmap(synfh, tclfh, qclfh, swap, k, r, trans, offsets);

    return synmap.search(gfffh);
}


//' remove links that disagree with the synteny map
//'
//' @param syn      synteny map file name
//' @param hit      int file name
//' @param swap     reverse direction of synteny map (e.g. swap query and target) 
//' @param k        match fuziness, integer
//' @param r        score decay rate, 0 means no context, high means more context
//' @param trans    score transform methods, single character
//' @param offsets  6-element integer vector of [01] offsets
// [[Rcpp::export]]
Rcpp::CharacterVector c_filter(
    std::string syn,
    std::string hit,
    bool swap,
    int k,
    double r,
    char trans,
    std::vector<int> offsets
)
{
    FILE* synfh = fopen(syn.c_str(), "r");
    FILE* hitfh = fopen(hit.c_str(), "r");

    Synmap synmap(synfh, nullptr, nullptr, swap, k, r, trans, offsets);

    return synmap.filter(hitfh);
}

//' trace intervals across genomes
//'
//' @param syn     synteny map file name
//' @param gff     GFF file name
//' @param swap    reverse direction of synteny map (e.g. swap query and target) 
//' @param offsets 6-element integer vector of [01] offsets
// [[Rcpp::export]]
Rcpp::DataFrame c_map(
    std::string syn,
    std::string gff,
    bool swap,
    std::vector<int> offsets
)
{
    FILE* synfh = fopen(syn.c_str(), "r");
    FILE* gfffh = fopen(gff.c_str(), "r");

    Synmap synmap(synfh, nullptr, nullptr, swap, 0, 0, 'i', offsets);

    return synmap.map(gfffh);
}

//' count overlaps
//'
//' @param syn      synteny map file name
//' @param gff      GFF file name
//' @param swap     reverse direction of synteny map (e.g. swap query and target) 
//' @param offsets  6-element integer vector of [01] offsets
// [[Rcpp::export]]
Rcpp::DataFrame c_count(
    std::string syn,
    std::string gff,
    bool swap,
    std::vector<int> offsets
)
{
    FILE* synfh = fopen(syn.c_str(), "r");
    FILE* gfffh = fopen(gff.c_str(), "r");

    Synmap synmap(synfh, nullptr, nullptr, swap, 0, 0, 'i', offsets);

    return synmap.count(gfffh);
}
