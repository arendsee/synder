#include <string>

#include "global.h"
#include "synmap.h"

//' print all blocks with contiguous set ids
//'
//' @param syn      synteny map file name
//' @param swap     reverse direction of synteny map (e.g. swap query and target) 
//' @param trans    score transform methods, single character
//' @param k       match fuziness, integer
//' @param r       score decay rate, 0 means no context, high means more context
//' @param offsets  4-element integer vector of [01] offsets (start/stop
//'                 offsets for the synteny maps and the GFF)
// [[Rcpp::export]]
Rcpp::DataFrame c_dump (
    std::string syn,
    bool swap,
    char trans,
    int k,
    double r,
    std::vector<int> offsets
)
{
    Synmap synmap(syn, "", "", swap, k, r, trans, offsets);

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
//' @param offsets  4-element integer vector of [01] offsets (start/stop
//'                 offsets for the synteny maps and the GFF)
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
    Synmap synmap(syn, tcl, qcl, swap, k, r, trans, offsets);

    return synmap.search(gff);
}


//' remove links that disagree with the synteny map
//'
//' @param syn      synteny map file name
//' @param hit      int file name
//' @param swap     reverse direction of synteny map (e.g. swap query and target) 
//' @param k        match fuziness, integer
//' @param r        score decay rate, 0 means no context, high means more context
//' @param trans    score transform methods, single character
//' @param offsets  4-element integer vector of [01] offsets (start/stop
//'                 offsets for the synteny maps and the GFF)
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
    Synmap synmap(syn, "", "", swap, k, r, trans, offsets);

    return synmap.filter(hit);
}

//' trace intervals across genomes
//'
//' @param syn     synteny map file name
//' @param gff     GFF file name
//' @param swap    reverse direction of synteny map (e.g. swap query and target) 
//' @param offsets  4-element integer vector of [01] offsets (start/stop
//'                 offsets for the synteny maps and the GFF)
// [[Rcpp::export]]
Rcpp::DataFrame c_map(
    std::string syn,
    std::string gff,
    bool swap,
    std::vector<int> offsets
)
{
    Synmap synmap(syn, "", "", swap, 0, 0, 'i', offsets);

    return synmap.map(gff);
}

//' count overlaps
//'
//' @param syn      synteny map file name
//' @param gff      GFF file name
//' @param swap     reverse direction of synteny map (e.g. swap query and target) 
//' @param offsets  4-element integer vector of [01] offsets (start/stop
//'                 offsets for the synteny maps and the GFF)
// [[Rcpp::export]]
Rcpp::DataFrame c_count(
    std::string syn,
    std::string gff,
    bool swap,
    std::vector<int> offsets
)
{
    Synmap synmap(syn, "", "", swap, 0, 0, 'i', offsets);

    return synmap.count(gff);
}
