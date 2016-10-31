#ifndef __GENOME_H__
#define __GENOME_H__

#include "global.h"
#include "contig.h"
#include "many_contiguous_sets.h"

#include <map>
#include <stack>
#include <Rcpp.h>

class Genome {
private:
    std::string name;
    std::map<std::string, Contig*> contig;
    std::stack<Block> pool;

    Contig* add_contig(std::string contig_name);

public:

    Genome(std::string name);
    ~Genome();

    Rcpp::DataFrame as_data_frame();

    /** get contig by name, die if no matches */
    Contig* get_contig(std::string contig_name);

    /** create a new Block, new contigs are created as needed */
    Block* add_block(
        std::string contig_name,
        long        start,
        long        stop,
        double      score,
        char        strand
    );

    void set_contig_lengths(FILE* clfile);

    size_t size() {
        return contig.size();
    }

    std::string get_name() {
        return name;
    }

    void dump_blocks();

    /** Link blocks by next and prev stop and next and prev start */
    void link_block_corners();

    /** Link Contig to first and last blocks */
    void set_contig_corners();

    /** Set a unique index for each set of overlapping sequences */
    void set_overlap_group(long& offset);

    void link_adjacent_blocks();

    void merge_overlaps();

    void refresh();

    void link_contiguous_blocks(long k, size_t& setid);

    void transfer_contiguous_sets(Genome*);

    void validate();

};

#endif
