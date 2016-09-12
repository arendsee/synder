#ifndef __GENOME_H__
#define __GENOME_H__

#include "global.h"
#include "contig.h"
#include "block.h"

#include <map>
#include <stack>

class Genome {
private:
    std::string name;
    std::map<std::string, Contig*> contig;
    std::stack<Block> pool;

    Contig* add_contig(std::string contig_name);

public:

    Genome(std::string name);
    ~Genome();

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

    void set_contig_lengths();

    size_t size() {
        return contig.size();
    }

    std::string get_name() {
        return name;
    }

    void print(
        bool print_blocks=true,     // recursively print all blocks
        bool print_backwards=false  // print blocks sorted by stop position
    );

    /** Link blocks by next and prev stop and next and prev start */
    void link_block_corners();

    /** Link Contig to first and last blocks */
    void set_contig_corners();

    /** Set a unique index for each set of overlapping sequences */
    void set_overlap_group();

    /** Link each node to its adjacent neighbors
     *
     * Link blocks to nearest non-overlapping up and downstream blocks
     *
     * For example, given these for blocks:
     *  |---a---|
     *            |--b--|
     *             |----c----|
     *                     |---d---|
     *                               |---e---|
     * a->adj := (NULL, b)
     * b->adj := (a, e)
     * c->adj := (a, e)
     * d->adj := (a, e)
     * e->adj := (d, NULL)
     */
    void link_adjacent_blocks();
    /** Do one half of the job */
    void link_adjacent_blocks_directed(Contig* con, Direction d);

    void merge_overlaps();

    void link_contiguous_blocks(long k);

    void validate();

};

#endif
