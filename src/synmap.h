#ifndef __SYNMAP_H__
#define __SYNMAP_H__

#include "global.h"
#include "bound.h"
#include "genome.h"
#include "linked_interval.h"
#include "feature.h"
#include "types.h"


#include <iterator>
#include <list>
#include <array>
#include <Rcpp.h>


/** A pair of syntenically linked Genome objects  */
class Synmap
{
private:
    Genome* genome[2] = { nullptr, nullptr };
    FILE*   synfile   = nullptr;
    FILE*   tclfile   = nullptr;
    FILE*   qclfile   = nullptr;
    int     swap      = 0;
    long    k         = 0;
    double  r         = 0.001;
    char    trans     = 'i';

    // The {{ is needed to workaround a bug in old g++ compilers
    std::array<int,4> offsets = {{0,1,0,0}};

    // utility function for loading GFF files
    std::vector<Feature> gff2features(FILE* fh);

    // loads synfile and calls the below functions in proper order
    void load_blocks();

    // wrappers for Genome functions
    void link_blocks();

    /** Checks invariants - dies if anything goes wrong */
    void validate();

public:
    Synmap(
        FILE*  synfile,
        FILE*  tclfile,
        FILE*  qclfile,
        bool   swap,
        int    k,
        double r,
        char   trans,
        std::vector<int> offsets
    );

    ~Synmap();

    Contig* get_contig(size_t gid, const char* contig_name);

    Rcpp::DataFrame as_data_frame();

    Rcpp::DataFrame count(FILE* intfile);

    Rcpp::DataFrame map(FILE* intfile);

    Rcpp::DataFrame search(FILE* intfile);

    Rcpp::CharacterVector filter(FILE* hitfile);

};

#endif
