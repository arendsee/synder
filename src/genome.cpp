#include "genome.h"

Genome::Genome(std::string new_name)
{
    name = new_name;
}

Genome::~Genome()
{
    for (auto &pair : contig)
    {
        delete pair.second;
    }
}

Contig* Genome::add_contig(std::string contig_name)
{
    Contig* con;
    auto it = contig.find(contig_name);
    if (it == contig.end())
    {
        con = new Contig(contig_name, this);
        contig[contig_name] = con;
    }
    else
    {
        con = (*it).second;
    }
    return con;
}

Contig* Genome::get_contig(std::string name)
{
    Contig* con;
    try
    {
        con = contig.at(name);
    }
    catch (const std::out_of_range& e)
    {
        con = NULL;
    }
    return con;
}

Block* Genome::add_block(
    std::string contig_name,
    long start,
    long stop,
    double score,
    char strand = '+'
)
{
    Contig* con = add_contig(contig_name);

    Block blk;

    blk.over   = NULL;
    blk.cor[0] = NULL;
    blk.cor[1] = NULL;
    blk.cor[2] = NULL;
    blk.cor[3] = NULL;
    blk.adj[0] = NULL;
    blk.adj[1] = NULL;
    blk.cor[0] = NULL;
    blk.cor[1] = NULL;
    blk.cnr[0] = NULL;
    blk.cnr[1] = NULL;
    blk.cset   = NULL;

    blk.parent = con;
    blk.pos[0] = start;
    blk.pos[1] = stop;
    blk.score  = score;
    blk.strand = strand;
    blk.linkid = pool.size() + 1;

    pool.push(blk);

    Block* blk_ptr = &pool.top();

    con->block.push_back(blk_ptr);

    return blk_ptr;
}

void Genome::set_contig_lengths()
{
    // TODO process file if given, to replace the default of 1000000000
    // Contig * con;
    // for (auto &pair : contig)
    // {
    //     con = &pair.second;
    // }
}

void Genome::print(bool print_blocks, bool forward)
{
    fprintf(stderr, ">%s size=%lu\n", name.c_str(), contig.size());
    for (auto &pair : contig)
    {
        pair.second->print(forward, print_blocks);
    }
}

void Genome::link_block_corners()
{
    size_t N;
    for (auto &pair : contig)
    {
        N = pair.second->block.size();

        Block** blocks = &pair.second->block[0];

        // sort by stop
        Contig::sort_blocks(blocks, N, true);
        for (size_t i = 0; i < N; i++)
        {
            blocks[i]->cor[2] = (i == 0)     ? NULL : blocks[i - 1];
            blocks[i]->cor[3] = (i == N - 1) ? NULL : blocks[i + 1];
        }
        // sort by start
        Contig::sort_blocks(blocks, N, false);
        for (size_t i = 0; i < N; i++)
        {
            blocks[i]->cor[0] = (i == 0)     ? NULL : blocks[i - 1];
            blocks[i]->cor[1] = (i == N - 1) ? NULL : blocks[i + 1];
        }
    }
}

void Genome::set_contig_corners()
{
    Contig * con;
    size_t k;
    for (auto &pair : contig)
    {
        con = pair.second;
        for (size_t i = 0; i < 4; i++)
        {
            k = i % 2 == 0 ? 0 : con->block.size() - 1;
            con->cor[i] = con->block[k];
            while (con->cor[i]->cor[i] != NULL)
            {
                con->cor[i] = con->cor[i]->cor[i];
            }
        }
    }
}


void Genome::set_overlap_group()
{

    // Holds current overlapping group id
    size_t grpid = 1;
    // Needed for determining overlaps and thus setids
    long maximum_stop = 0;
    // The stop position of the current interval
    long this_stop = 0;
    // Current Block in linked list
    Block* blk;
    Contig* con;

    // Loop through each contig in the query genome
    // i := contig id
    for (auto &pair : contig)
    {
        con = pair.second;
        maximum_stop = 0;
        // Loop through each Block in the linked list
        blk = con->cor[0];
        for (; blk != NULL; blk = blk->cor[1])
        {
            this_stop = blk->pos[1];
            // If the start is greater than the maximum stop, then the block is in
            // a new adjacency group. For this to work, Contig->block must be
            // sorted by start. This sort is performed in build_tree.
            if (blk->pos[0] > maximum_stop)
            {
                grpid++;
            }
            if (this_stop > maximum_stop)
            {
                maximum_stop = this_stop;
            }
            blk->grpid = grpid;
        }
        // increment to break adjacency between contigs and genomes
        grpid++;
    }
}


void Genome::link_adjacent_blocks_directed(Contig* con, Direction d)
{
    // In diagrams:
    // <--- indicates a hi block
    // ---> indicates a lo block
    // All diagrams and comments relative to the d==HI direction

    if (con->cor[0] == NULL || con->cor[1] == NULL || con->cor[2] == NULL || con->cor[3] == NULL)
    {
        fprintf(stderr, "Contig head must be set before link_adjacent_blocks is called\n");
        fprintf(stderr, "genome=(%s) contig=(%s)\n", con->parent->name.c_str(), con->name.c_str());
        exit(EXIT_FAILURE);
    }

    // Transformed indices for Block->cor and Contig->cor
    int idx_a = (!d * 2) + !d; // - 0 previous/first element by start
    int idx_b = (!d * 2) +  d; // - 1 next/last element by start
    int idx_c = ( d * 2) + !d; // - 2 previous/first element by stop
    int idx_d = ( d * 2) +  d; // - 3 next/last element by stop

    Block *lo, *hi;

    lo = con->cor[idx_c]; // first element by stop
    hi = con->cor[idx_a]; // first element by start

    while (hi != NULL)
    {

        //       --->
        // <---
        // OR
        // --->
        //   <---
        // This should should occur only at the beginning
        if (REL_LE(hi->pos[!d], lo->pos[d], d))
        {
            hi->adj[!d] = NULL;
            hi = hi->cor[idx_b];
        }

        //  lo     next
        // ---->  ---->
        //               <---
        // If next is closer, and not overlapping the hi, increment lo
        // You increment lo until it is adjacent to the current hi
        else if (REL_LT(lo->cor[idx_d]->pos[d], hi->pos[!d], d))
        {
            lo = lo->cor[idx_d];
        }

        // --->
        //      <---
        // The current lo is next to, and not overlapping, current hi
        else
        {
            hi->adj[!d] = lo;
            hi = hi->cor[idx_b];
        }
    }
}

void Genome::link_adjacent_blocks()
{
    Contig * con;
    for (auto &pair : contig)
    {
        con = pair.second;
        link_adjacent_blocks_directed(con, HI);
        link_adjacent_blocks_directed(con, LO);
    }
}

void Genome::merge_overlaps()
{
    for (auto &pair : contig)
    {
        pair.second->merge_doubly_overlapping_blocks();
    }
}

void Genome::link_contiguous_blocks(long k)
{

    Block * blk;

    std::list<ContiguousSet*> csets;
    std::list<ContiguousSet*>::iterator iter;
    ContiguousSet * cset;

    for (auto &pair : contig)
    {

        // Initialize the first block in the scaffold
        blk = pair.second->cor[0];

        // Initialize first ContiguousSet
        csets.clear();
        csets.push_front(init_ContiguousSet(blk));

        for (blk = blk->cor[1]; blk != NULL; blk = blk->cor[1])
        {
            for (auto &cset : csets)
            {
                // if block has joined a set
                if (add_block_to_ContiguousSet(cset, blk, k))
                    goto added;

                // if set terminates
                if (strictly_forbidden(cset->ends[1], blk, k))
                    iter = csets.erase(iter);
            }
            // if block fits in no set, create a new one
            csets.push_front(init_ContiguousSet(blk));
added:
            {}
        }
    }

    size_t setid = 0;
    for (auto &pair : contig)
    {
        cset = pair.second->cset;
        // rewind - TODO - build the csets so this isn't necessary
        while (cset->prev != NULL)
        {
            cset = cset->prev;
        }
        pair.second->cset = cset;
        // TODO - remove setids. I don't use them for anything but debugging
        for (; cset != NULL; cset = cset->next)
        {
            setid++;
            cset->id = setid;
            cset->over->id = setid;
        }
    }
}

void Genome::validate()
{
#define ASSERT_CON(t, con)                                   \
        if(!(t)){                                            \
          is_good=false;                                     \
          fprintf(stderr, "Assert failed: `" #t "`\n");      \
          if(blk != NULL)                                    \
            con->print(false, false);                        \
        }
#define ASSERT_BLK(t, blk)                                   \
        if(!(t)){                                            \
          is_good=false;                                     \
          fprintf(stderr, "Assert failed: `" #t "`\n");      \
          if(blk != NULL){print_Block(blk);}                 \
          else {fprintf(stderr, "NULL\n");}                  \
        }

    size_t  nblks   = 0;
    Contig* con     = NULL;
    Block*  blk     = NULL;
    bool    is_good = true;

    for (auto &pair : contig)
    {
        con = pair.second;

        ASSERT_CON(con->cor[0] != NULL, con);
        ASSERT_CON(con->cor[1] != NULL, con);
        ASSERT_CON(con->cor[2] != NULL, con);
        ASSERT_CON(con->cor[3] != NULL, con);

        blk = con->cor[0];
        nblks = 0;
        for (; blk != NULL; blk = blk->cor[1])
        {
            if (!(blk->pos[1] < con->length))
            {
                fprintf(stderr,
                        "WARNING: stop greater than contig length: %zu vs %zu\n",
                        blk->pos[1], con->length);
            }

            nblks++;

            ASSERT_BLK(blk->over != NULL, blk);
            ASSERT_BLK(blk->cset != NULL, blk);
            ASSERT_BLK(blk == blk->over->over, blk);
            ASSERT_BLK(blk->cset->id == blk->over->cset->id, blk);
            ASSERT_BLK(blk->cset->over == blk->over->cset, blk);
            ASSERT_BLK(blk->score == blk->over->score, blk);

            ASSERT_BLK(blk->pos[0] >= con->cor[0]->pos[0], blk);
            ASSERT_BLK(blk->pos[0] <= con->cor[1]->pos[0], blk);
            ASSERT_BLK(blk->pos[1] >= con->cor[2]->pos[1], blk);
            ASSERT_BLK(blk->pos[1] <= con->cor[3]->pos[1], blk);

            if (blk->cor[0] != NULL)
            {
                ASSERT_BLK(blk->pos[0] >= blk->cor[0]->pos[0], blk);
                ASSERT_BLK(blk->cor[0]->cor[1]->pos[0] == blk->pos[0], blk);
            }
            if (blk->cor[1] != NULL)
            {
                ASSERT_BLK(blk->pos[0] <= blk->cor[1]->pos[0], blk);
                ASSERT_BLK(blk->cor[1]->cor[0]->pos[0] == blk->pos[0], blk);
            }
            if (blk->cor[2] != NULL)
            {
                ASSERT_BLK(blk->pos[1] >= blk->cor[2]->pos[1], blk);
                ASSERT_BLK(blk->cor[2]->cor[3]->pos[1] == blk->pos[1], blk);
            }
            if (blk->cor[3] != NULL)
            {
                ASSERT_BLK(blk->pos[1] <= blk->cor[3]->pos[1], blk);
                ASSERT_BLK(blk->cor[3]->cor[2]->pos[1] == blk->pos[1], blk);
            }

            // grpid == 0 only if unset
            ASSERT_BLK(blk->grpid != 0, blk);

            for (size_t i = 0; i < 2; i++)
            {
                if (blk->cnr[i] != NULL)
                {
                    ASSERT_BLK(blk->grpid != blk->cnr[i]->grpid, blk);
                    ASSERT_BLK(blk->cset == blk->cnr[i]->cset, blk);
                    ASSERT_BLK(blk->cnr[i]->over->cnr[!i]->over == blk, blk);
                }
            }

        }
        // blocks may be deleted, so nblks == con->size may not hold
        // however no blocks should ever be added
        ASSERT_CON(nblks <= con->block.size(), con);
    }

    if (! is_good)
        exit(EXIT_FAILURE);

#undef ASSERT_BLK
#undef ASSERT_CON
}
