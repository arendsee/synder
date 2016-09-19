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

Block* Genome::add_block(std::string contig_name, long start, long stop, double score, char strand = '+')
{
    Contig* con = add_contig(contig_name);

    Block blk(start, stop, score, strand, con, pool.size() + 1);

    pool.push(blk);

    Block* blk_ptr = &pool.top();

    con->block->add(blk_ptr);

    return blk_ptr;
}

void Genome::set_contig_lengths(FILE* clfile)
{
    if(clfile != NULL){
        // read loop variables
        int     status = 0;
        char*   line   = NULL;
        size_t  len    = 0;
        ssize_t read;

        char contig_name[NAME_BUFFER_SIZE];
        long contig_length;

        Contig* con;
        while ((read = getline(&line, &len, clfile)) != EOF)
        {
            status = sscanf(line, "%s %ld", contig_name, &contig_length);
            if(status == EOF){
                fprintf(stderr, "Failed to read contig length file\n");
                exit(EXIT_FAILURE);
            }
            con = get_contig(contig_name);
            if(con != NULL){
                con->length = contig_length;
            }
        }
        free(line);
    }
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
    Contig* con;
    for (auto &pair : contig)
    {
        con = pair.second;
        con->block->link_block_corners();
    }
}

void Genome::set_contig_corners()
{
    Contig* con;
    for (auto &pair : contig)
    {
        con = pair.second;
        con->set_contig_corners();
    }
}


void Genome::set_overlap_group()
{

    Contig* con;

    // Loop through each contig in the query genome
    // i := contig id
    for (auto &pair : contig)
    {
        con = pair.second;
        con->block->set_overlap_group();
    }
}

void Genome::link_adjacent_blocks()
{
    Contig* con;
    for (auto &pair : contig)
    {
        con = pair.second;
        con->block->link_adjacent_blocks();
    }
}

void Genome::merge_overlaps()
{
    for (auto &pair : contig)
    {
        pair.second->block->merge_overlaps();
    }
}

void Genome::link_contiguous_blocks(long k)
{
    size_t setid = 0;
    for (auto &pair : contig)
    {
        pair.second->link_contiguous_blocks(k, setid);
    }
}

void Genome::validate()
{
/*
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
*/
}
