#include "genome.h"

Genome::Genome(std::string new_name)
{
    name = new_name;
}

Genome::~Genome()
{
    for (auto &pair : contig) {
        delete pair.second;
    }
}

Contig* Genome::add_contig(std::string contig_name)
{
    Contig* con;
    auto it = contig.find(contig_name.c_str());
    if (it == contig.end()) {
        con = new Contig(name.c_str(), contig_name.c_str());
        contig[contig_name.c_str()] = con;
    } else {
        con = (*it).second;
    }
    return con;
}

Contig* Genome::get_contig(std::string t_name)
{
    Contig* con;
    try {
        con = contig.at(t_name);
    } catch (const std::out_of_range& e) {
        cerr << "Index error in " << __func__ << endl;
        con = nullptr;
    }
    return con;
}

Block* Genome::add_block(std::string contig_name, long start, long stop, double score, char strand = '+')
{
    Contig* con = add_contig(contig_name.c_str());

    pool.push(
        Block(start, stop, score, strand, &con->feat, pool.size() + 1)
    );

    Block* blk_ptr = &pool.top();

    con->block.add(blk_ptr);

    return blk_ptr;
}

void Genome::set_contig_lengths(FILE* clfile)
{
    if(clfile != nullptr) {
        // read loop variables
        int     status = 0;
        char*   line   = nullptr;
        size_t  len    = 0;
        ssize_t read;

        char contig_name[NAME_BUFFER_SIZE];
        long contig_length;

        Contig* con;
        while ((read = getline(&line, &len, clfile)) != EOF) {
            status = sscanf(line, "%s %ld", contig_name, &contig_length);
            if(status == EOF) {
                fprintf(stderr, "Failed to read contig length file\n");
                exit(EXIT_FAILURE);
            }
            con = get_contig(contig_name);
            if(con != nullptr) {
                con->set_length(contig_length);
            }
        }
        free(line);
    }
}

void Genome::dump_blocks()
{
    for (auto &pair : contig) {
        pair.second->block.print();
    }
}

void Genome::print(bool print_blocks, bool forward)
{
    fprintf(stderr, ">%s size=%lu\n", name.c_str(), contig.size());
    for (auto &pair : contig) {
        pair.second->print(forward, print_blocks);
    }
}

void Genome::link_block_corners()
{
    for (auto &pair : contig) {
        pair.second->block.link_block_corners();
    }
}

void Genome::set_contig_corners()
{
    for (auto &pair : contig) {
        pair.second->block.link_corners();
    }
}

void Genome::set_overlap_group()
{
    for (auto &pair : contig) {
        pair.second->block.set_overlap_group();
    }
}

void Genome::link_adjacent_blocks()
{
    for (auto &pair : contig) {
        pair.second->block.link_adjacent_blocks();
    }
}

void Genome::merge_overlaps()
{
    for (auto &pair : contig) {
        pair.second->block.merge_overlaps();
    }
}

void Genome::link_contiguous_blocks(long k, size_t& setid)
{
    Contig * con;
    Block * first_blk;
    for (auto &pair : contig) {
        con = pair.second;
        first_blk = con->block.front();
        con->cset.link_contiguous_blocks(first_blk, k, setid);
    }
}

void Genome::validate()
{
    #define ASSERT_CON(t)                                        \
            if(!(t)){                                            \
              is_good=false;                                     \
              fprintf(stderr, "Assert failed: `" #t "`\n");      \
            }
    #define ASSERT_BLK(t)                                        \
            if(!(t)){                                            \
              is_good=false;                                     \
              fprintf(stderr, "Assert failed: `" #t "`\n");      \
            }

        bool is_good = true;

        for (auto &pair : contig)
        {
            Contig* con = pair.second;

            ASSERT_CON(con->block.corner(0) != nullptr);
            ASSERT_CON(con->block.corner(1) != nullptr);
            ASSERT_CON(con->block.corner(2) != nullptr);
            ASSERT_CON(con->block.corner(3) != nullptr);

            Block* blk = con->block.corner(0);
            for (; blk != nullptr; blk = blk->corner(1))
            {
                ASSERT_BLK(blk->stop() <= con->feat.parent_length)

                ASSERT_BLK(blk == blk->over->over);
                ASSERT_BLK(blk->cset->id   == blk->over->cset->id);
                ASSERT_BLK(blk->cset->over == blk->over->cset);
                ASSERT_BLK(blk->score      == blk->over->score);

                ASSERT_BLK(blk->pos[0] >= con->block.corner(0)->pos[0]);
                ASSERT_BLK(blk->pos[1] >= con->block.corner(2)->pos[1]);

                // ASSERT_BLK(blk->pos[0] <= con->block.corner(1)->pos[0]);
                // ASSERT_BLK(blk->pos[1] <= con->block.corner(3)->pos[1]);

                if (blk->corner(0) != nullptr)
                {
                    ASSERT_BLK(blk->pos[0] >= blk->corner(0)->pos[0]);
                    ASSERT_BLK(blk->corner(0)->corner(1)->pos[0] == blk->pos[0]);
                }
                if (blk->corner(1) != nullptr)
                {
                    ASSERT_BLK(blk->pos[0] <= blk->corner(1)->pos[0]);
                    ASSERT_BLK(blk->corner(1)->corner(0)->pos[0] == blk->pos[0]);
                }
                if (blk->corner(2) != nullptr)
                {
                    ASSERT_BLK(blk->pos[1] >= blk->corner(2)->pos[1]);
                    ASSERT_BLK(blk->corner(2)->corner(3)->pos[1] == blk->pos[1]);
                }
                if (blk->corner(3) != nullptr)
                {
                    ASSERT_BLK(blk->pos[1] <= blk->corner(3)->pos[1]);
                    ASSERT_BLK(blk->corner(3)->corner(2)->pos[1] == blk->pos[1]);
                }

                // grpid == 0 only if unset
                ASSERT_BLK(blk->grpid != 0);

                for (size_t i = 0; i < 2; i++)
                {
                    if (blk->cnr[i] != nullptr)
                    {
                        ASSERT_BLK(blk->grpid != blk->cnr[i]->grpid);
                        ASSERT_BLK(blk->cset == blk->cnr[i]->cset);
                        ASSERT_BLK(blk->cnr[i]->over->cnr[!i]->over == blk);
                    }
                }

            }
        }


        if (! is_good){
            throw "Synmap invariants did";
        }

    #undef ASSERT_BLK
    #undef ASSERT_CON
}
