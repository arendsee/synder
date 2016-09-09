#include "synmap.h"

Synmap::Synmap(FILE * synfile, int swap, long k, char trans, bool validate_on)
{
    if (synfile == NULL)
    {
        fprintf(stderr, "NULL synteny input file (%s:%d in %s)", __FILE__, __LINE__, __func__);
        exit(EXIT_FAILURE);
    }


    // for swap == 0, the query genome is the first in the synmap
    // for swap == 1, it is the second
    int query  = swap;
    int target = !swap;

    // read loop variables
    size_t  line_no         = 0;
    size_t  unloaded_blocks = 0;
    int     status          = 0;
    char    *line           = NULL;
    size_t  len             = 0;
    ssize_t read;

    // Contig name
    char    seqid[128];
    // contig id
    size_t  conid;
    // length of contig in bases (i.e. scaffold or chromosome length)
    size_t  contig_length;
    // number of blocks in contig
    size_t  nblocks;
    // number of contigs in current genome
    size_t  ncontigs;

    // NOTE: The synteny database files are always 0-based. So subtracting
    // global_in_base from the intervals is not necessary.

    // A dummy variable searched for at the end of each sscanf format sequence.
    // It will be matched only if there are too many arguments.
    char dummy;
    for (int i = 0; i < 2; i++)
    {
        while ((read = getline(&line, &len, synfile)) != EOF)
        {
            line_no++;
            size_t loc = i == query ? 0 : 1;
            if (line[0] == '>')
            {
                status = sscanf(line, "> %s %zu %c", seqid, &ncontigs, &dummy);
                check_args(line_no, status, 2);
                genome[loc] = init_Genome(seqid, ncontigs);
            }
            else if (line[0] == '@')
            {
                break;
            }
            else if (line[0] == '$')
            {
                status = sscanf(line, "$ %zu %zu %s %zu %c\n", &conid, &nblocks, seqid, &contig_length, &dummy);
                check_args(line_no, status, 4);
                genome[loc]->contig[conid] = init_Contig(seqid, nblocks, contig_length, genome[loc]);
                unloaded_blocks += nblocks;
            }
            else
            {
                fprintf(stderr, "Incorrect file format, line %zu\n", line_no);
                fprintf(stderr, "Offending line:\n%s", line);
                exit(EXIT_FAILURE);
            }
        }
    }

    size_t qcon_id, qblk_id, tcon_id, tblk_id;
    long qstart, qstop, tstart, tstop;
    double score;
    char strand;

    Block *qblk, *tblk;
    Contig *tcon, *qcon;

    for (
            line_no = 0;
            read != EOF;
            line_no++, read = getline(&line, &len, synfile)
        )
    {
        if (line[0] != '$')
            continue;

        unloaded_blocks -= 2;
        status = sscanf(line, "$ %zu %zu %zu %zu %zu %zu %zu %zu %lf %c %c\n",
                        &qcon_id, &qblk_id, &qstart, &qstop,
                        &tcon_id, &tblk_id, &tstart, &tstop, &score, &strand, &dummy);
        check_args(line_no, status, 10);

        qcon = genome[query]->contig[qcon_id];
        tcon = genome[target]->contig[tcon_id];

        if (qstart > qstop || tstart > tstop)
        {
            fprintf(stderr, "start must be less than stop on line %zu\n", line_no);
            fprintf(stderr, "offending line:\n%s\n", line);
            exit(EXIT_FAILURE);
        }
        if (tcon->length <= tstop)
        {
            fprintf(stderr, "stop must be less than contig length on line %zu\n", line_no);
            fprintf(stderr, "conid=%zu blkid=%zu pos=(%zu, %zu) conlen=%zu\n",
                    tcon_id, tblk_id, tstart, tstop, tcon->length);
            fprintf(stderr, "offending line:\n%s\n", line);
            exit(EXIT_FAILURE);
        }
        // don't exceed the specified number of Contig in Genome
        if (qcon_id >= genome[query]->size
                || tcon_id >= genome[target]->size)
        {
            fprintf(stderr, "too few contigs specified\n");
            exit(EXIT_FAILURE);
        }
        // don't exceed the specified size of the Contig block arrays
        if (qblk_id >= qcon->size ||
                tblk_id >= tcon->size)
        {
            fprintf(stderr, "too few blocks specified\n");
            exit(EXIT_FAILURE);
        }

        qblk = &qcon->block[qblk_id];
        tblk = &tcon->block[tblk_id];

        switch (trans)
        {
        case 'l':
            score = -1 * log(score);
            break;
        case 'd':
            score = score * MIN((tstop - tstart + 1), (qstop - qstart + 1));
            break;
        case 'p':
            score = score * MIN((tstop - tstart + 1), (qstop - qstart + 1)) / 100.0;
            break;
        case 'i':
            // no transformation
            break;
        default:
            fprintf(stderr, "Unexpected transformation '%c'\n", trans);
            exit(EXIT_FAILURE);
            break;
        }

        set_Block(qblk, qstart, qstop, score, '+',    qcon, tblk, line_no);
        set_Block(tblk, tstart, tstop, score, strand, tcon, qblk, line_no);

    }
    free(line);

    /* each block specified has been loaded */
    if (unloaded_blocks != 0)
    {
        fprintf(stderr, "Unequal number of specified and actual blocks\n");
        exit(EXIT_FAILURE);
    }

    // The following must be run in order
    link_block_corners();
    set_contig_corners();
    merge_all_doubly_overlapping_blocks();
    set_overlap_group();
    link_adjacent_blocks();
    link_contiguous_blocks(k);
    validate();
}

Contig * Synmap::get_contig(size_t gid, size_t cid){
    return genome[gid]->contig[cid];
}

Genome * Synmap::get_genome(size_t gid){
    return genome[gid];
}

void Synmap::check_args(size_t line_no, size_t nargs, size_t correct_nargs)
{
    if (nargs < correct_nargs)
    {
        fprintf(stderr, "Either too few fields on input line %zu,", line_no);
        fprintf(stderr, "or they are of the wrong type\n");
        fprintf(stderr, "Found %zu arguments, expected %zu.\n", nargs, correct_nargs);
        exit(EXIT_FAILURE);
    }
    else if (nargs > correct_nargs)
    {
        fprintf(stderr, "Too many fields on input line %zu\n", line_no);
        exit(EXIT_FAILURE);
    }
}

Synmap::~Synmap()
{
    // TODO
    free_Genome(genome[0]);
    free_Genome(genome[1]);
}

void Synmap::print(bool forward)
{
    // only print the query Genome, the print_verbose_Block function will print
    // the target information as well
    fprintf(
        stderr,
        "--- Query=(%s, %zu), Target=(%s, %zu)\n",
        genome[0]->name,
        genome[0]->size,
        genome[1]->name,
        genome[1]->size
    );
    fprintf(stderr, "---------------------------------------------------------\n");
    fprintf(stderr, "Target contigs:\n");
    print_Genome(genome[1], forward, false);
    fprintf(stderr, "---------------------------------------------------------\n");
    fprintf(stderr, "Query contigs and blocks:\n");
    print_Genome(genome[0], forward, true);
}

void Synmap::dump_blocks()
{
    for (size_t i = 0; i < genome[0]->size; i++)
    {
        Block * blk = get_contig(0, i)->cor[0];
        for (; blk != NULL; blk = blk->cor[1])
        {
            print_Block(blk);
        }
    }
}

void Synmap::link_block_corners()
{
    Contig* con;
    size_t N;
    for (size_t g = 0; g <= 1; g++)
    {
        for (size_t c = 0; c < get_genome(g)->size; c++)
        {
            con = get_contig(g, c);
            N = con->size;

            Block** blocks = (Block**)malloc(N * sizeof(Block*));
            for (size_t i = 0; i < N; i++)
            {
                blocks[i] = &con->block[i];
            }

            // sort by stop
            sort_blocks(blocks, N, true);
            for (size_t i = 0; i < N; i++)
            {
                blocks[i]->cor[2] = (i == 0)     ? NULL : blocks[i - 1];
                blocks[i]->cor[3] = (i == N - 1) ? NULL : blocks[i + 1];
            }
            // sort by start
            sort_blocks(blocks, N, false);
            for (size_t i = 0; i < N; i++)
            {
                blocks[i]->cor[0] = (i == 0)     ? NULL : blocks[i - 1];
                blocks[i]->cor[1] = (i == N - 1) ? NULL : blocks[i + 1];
            }

            free(blocks);
        }
    }
}

void Synmap::set_contig_corners()
{
    Contig* con;
    int k;
    for (size_t g = 0; g < 2; g++)
    {
        for (size_t c = 0; c < get_genome(g)->size; c++)
        {
            con = get_contig(g, c);
            for (size_t i = 0; i < 4; i++)
            {
                k = i % 2 == 0 ? 0 : con->size - 1;
                con->cor[i] = &con->block[k];
                while (con->cor[i]->cor[i] != NULL)
                {
                    con->cor[i] = con->cor[i]->cor[i];
                }
            }
        }
    }
}

void Synmap::set_overlap_group()
{

    // Holds current overlapping group id
    size_t grpid = 1;
    // Needed for determining overlaps and thus setids
    long maximum_stop = 0;
    // The stop position of the current interval
    long this_stop = 0;
    // Current Block in linked list
    Block* blk;

    // Loop through target and query genomes
    // g := genome id (0 is query, 1 is target)
    for (size_t g = 0; g <= 1; g++)
    {
        // Loop through each contig in the query genome
        // i := contig id
        for (size_t i = 0; i < get_genome(g)->size; i++)
        {
            maximum_stop = 0;
            // Loop through each Block in the linked list
            blk = get_contig(g, i)->cor[0];
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
}


void Synmap::link_adjacent_blocks_directed(Contig* con, Direction d)
{
    // In diagrams:
    // <--- indicates a hi block
    // ---> indicates a lo block
    // All diagrams and comments relative to the d==HI direction

    if (con->cor[0] == NULL || con->cor[1] == NULL || con->cor[2] == NULL || con->cor[3] == NULL)
    {
        fprintf(stderr, "Contig head must be set before link_adjacent_blocks is called\n");
        fprintf(stderr, "genome=(%s) contig=(%s)\n", con->parent->name, con->name);
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
void Synmap::link_adjacent_blocks()
{
    for (size_t genid = 0; genid <= 1; genid++)
    {
        for (size_t conid = 0; conid < get_genome(genid)->size; conid++)
        {
            link_adjacent_blocks_directed(get_contig(genid, conid), HI);
            link_adjacent_blocks_directed(get_contig(genid, conid), LO);
        }
    }
}

void Synmap::merge_all_doubly_overlapping_blocks()
{
    Contig * con;
    for (size_t i = 0; i < get_genome(0)->size; i++)
    {
        con = get_contig(0, i);
        merge_doubly_overlapping_blocks(con);
    }
}


void Synmap::link_contiguous_blocks(long k)
{

    Contig * con;
    Block * blk;

    std::list<ContiguousSet*> csets;
    std::list<ContiguousSet*>::iterator iter;
    ContiguousSet * cset;

    for (size_t i = 0; i < get_genome(0)->size; i++)
    {

        // Initialize the first block in the scaffold
        blk = get_contig(0, i)->cor[0];

        // Initialize first ContiguousSet
        csets.clear();
        csets.push_front(init_ContiguousSet(blk));

        for (blk = blk->cor[1]; blk != NULL; blk = blk->cor[1])
        {
            for(iter = csets.begin(); iter != csets.end(); ++iter)
            {
                cset = *iter;

                // if block has joined a set
                if (add_block_to_ContiguousSet(cset, blk, k))
                    goto added;

                // if set terminates
                if (strictly_forbidden(cset->ends[1], blk, k))
                    iter = csets.erase(iter);
            }
                // if block fits in no set, create a new one
                csets.push_front(init_ContiguousSet(blk));
            added: {}
        }
    }

    size_t setid = 0;
    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 0; j < get_genome(i)->size; j++)
        {
            con = get_contig(i, j);
            // rewind - TODO - build the csets so this isn't necessary
            while (con->cset->prev != NULL)
            {
                con->cset = con->cset->prev;
            }
            if (i == 0)
            {
                // TODO - remove setids. I don't use them for anything but debugging
                for (cset = con->cset; cset != NULL; cset = cset->next)
                {
                    setid++;
                    cset->id = setid;
                    cset->over->id = setid;
                }
            }
        }
    }
}

void Synmap::validate()
{
#define ASSERT_CON(t, con)                                   \
        if(!(t)){                                            \
          is_good=false;                                     \
          fprintf(stderr, "Assert failed: `" #t "`\n");      \
          if(blk != NULL)                                    \
            print_Contig(con, false, false);                 \
        }
#define ASSERT_BLK(t, blk)                                   \
        if(!(t)){                                            \
          is_good=false;                                     \
          fprintf(stderr, "Assert failed: `" #t "`\n");      \
          if(blk != NULL){print_Block(blk);}                 \
          else {fprintf(stderr, "NULL\n");}                  \
        }
    size_t  gid     = 0;
    size_t  cid     = 0;
    size_t  nblks   = 0;
    Contig* con     = NULL;
    Block*  blk     = NULL;
    bool    is_good = true;

    for (gid = 0; gid < 2; gid++)
    {
        for (cid = 0; cid < get_genome(gid)->size; cid++)
        {
            con = get_contig(gid, cid);

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
            ASSERT_CON(nblks <= con->size, con);
        }
    }
    if (! is_good)
    {
        print(true);
        exit(EXIT_FAILURE);
    }
#undef ASSERT_BLK
#undef ASSERT_CON
}

void Synmap::count(FILE * intfile)
{
  char seqname[128];
  size_t count;
  size_t chrid;
  long start, stop;
  while ((fscanf(intfile,
                 "%zu %*s %*s %li %li %*s %*c %*s %s\n",
                 &chrid, &start, &stop, seqname)) != EOF)
  {
    check_in_offset(start, stop);
    start -= Offsets::in_start;
    stop  -= Offsets::in_stop;

    count = count_overlaps(get_contig(0, chrid), start, stop);
    printf("%s\t%zu\n", seqname, count);
  }
}

void Synmap::map(FILE * intfile)
{
  char seqname[128];
  size_t chrid;
  long start, stop;
  ResultContig * rc;
  Block *qblk, *tblk;
  bool missing;
  while ((fscanf(intfile,
                 "%zu %*s %*s %zu %zu %*s %*c %*s %s\n",
                 &chrid, &start, &stop, seqname)) != EOF)
  {
    check_in_offset(start, stop);
    start -= Offsets::in_start;
    stop  -= Offsets::in_stop;

    rc = get_region(get_contig(0, chrid), start, stop, false);
    missing = rc->inbetween || rc->leftmost || rc->rightmost;

    for (size_t i = 0; i < rc->size; i++) {
      qblk = rc->block[i];
      if (qblk != NULL) {
        tblk = qblk->over;
        printf("%s %s %zu %zu %d\n",
               seqname,
               tblk->parent->name,
               tblk->pos[0] + Offsets::out_start,
               tblk->pos[1] + Offsets::out_stop,
               missing
        );
      }
    }

    free(rc);
  }
}
