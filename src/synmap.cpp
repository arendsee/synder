#include "synmap.h"

Synmap::Synmap(Arguments& args)
{
    synfile = args.synfile;
    tclfile = args.tclfile;
    qclfile = args.qclfile;
    swap    = args.swap;
    k       = args.k;
    trans   = args.trans;

    load_blocks();
    set_contig_lengths();
    validate();
}

void Synmap::load_blocks()
{
    if (synfile == NULL) {
        fprintf(stderr, "NULL synteny input file (%s:%d in %s)", __FILE__, __LINE__, __func__);
        exit(EXIT_FAILURE);
    }

    genome[0] = new Genome("Q");
    genome[1] = new Genome("T");

    // read loop variables
    int     status = 0;
    char*   line   = NULL;
    size_t  len    = 0;
    ssize_t read;

    // Contig name
    char qseqid[NAME_BUFFER_SIZE];
    char tseqid[NAME_BUFFER_SIZE];
    char* seqid[2] = {qseqid, tseqid};
    double score;
    char strand;
    long start[2], stop[2];

    Block *qblk, *tblk;

    size_t i = swap ? 1 : 0;
    size_t j = swap ? 0 : 1;

    while ((read = getline(&line, &len, synfile)) != EOF) {

        // skip comments
        if (line[0] == '#')
            continue;

        status = sscanf(line, "%s %ld %ld %s %ld %ld %lf %c",
                        seqid[i], &start[i], &stop[i],
                        seqid[j], &start[j], &stop[j],
                        &score, &strand);

        if(status != 8) {
            fprintf(stderr, "Failed to read input line:\n%s", line);
            exit(EXIT_FAILURE);
        }

        switch (trans) {
            case 'l':
                score = -1 * log(score);
                break;
            case 'd':
                score = score * MIN((stop[0] - start[0] + 1), (stop[1] - start[1] + 1));
                break;
            case 'p':
                score = score * MIN((stop[0] - start[0] + 1), (stop[0] - start[0] + 1)) / 100.0;
                break;
            case 'i':
                // no transformation
                break;
            default:
                fprintf(stderr, "Unexpected transformation '%c'\n", trans);
                exit(EXIT_FAILURE);
                break;
        }

        qblk = genome[0]->add_block(seqid[0], start[0], stop[0], score, '+');
        tblk = genome[1]->add_block(seqid[1], start[1], stop[1], score, strand);

        // link homologs
        qblk->over = tblk;
        tblk->over = qblk;

    }
    free(line);

    // The following must be run in order
    link_block_corners();
    set_contig_corners();
    merge_doubly_overlapping_blocks();
    set_overlap_group();
    link_adjacent_blocks();
    link_contiguous_blocks(k);
}

Contig* Synmap::get_contig(size_t gid, char* contig_name)
{
    if (gid == 0 || gid == 1) {
        return genome[gid]->get_contig(contig_name);
    } else {
        return NULL;
    }
}

Synmap::~Synmap()
{
    delete genome[0];
    delete genome[1];
}

void Synmap::print(bool forward)
{
    // only print the query Genome, the print_verbose_Block function will print
    // the target information as well
    fprintf(
        stderr,
        "--- Query=(%s, %zu), Target=(%s, %zu)\n",
        genome[0]->get_name().c_str(),
        genome[0]->size(),
        genome[1]->get_name().c_str(),
        genome[1]->size()
    );
    fprintf(stderr, "---------------------------------------------------------\n");
    fprintf(stderr, "Target contigs:\n");
    genome[1]->print(forward, false);
    fprintf(stderr, "---------------------------------------------------------\n");
    fprintf(stderr, "Query contigs and blocks:\n");
    genome[0]->print(forward, true);
}

void Synmap::dump_blocks()
{
    // TODO implement
}

void Synmap::set_contig_lengths()
{
    genome[0]->set_contig_lengths(qclfile);
    genome[1]->set_contig_lengths(tclfile);
}

void Synmap::link_block_corners()
{
    genome[0]->link_block_corners();
    genome[1]->link_block_corners();
}

void Synmap::set_contig_corners()
{
    genome[0]->set_contig_corners();
    genome[1]->set_contig_corners();
}

void Synmap::set_overlap_group()
{
    genome[0]->set_overlap_group();
    genome[1]->set_overlap_group();
}

void Synmap::link_adjacent_blocks()
{
    genome[0]->link_adjacent_blocks();
    genome[1]->link_adjacent_blocks();
}

void Synmap::merge_doubly_overlapping_blocks()
{
    genome[0]->merge_overlaps();
}

void Synmap::link_contiguous_blocks(long k)
{
    genome[0]->link_contiguous_blocks(k);
}

void Synmap::validate()
{
    genome[0]->validate();
    genome[1]->validate();
}


bool Synmap::process_gff(FILE* intfile, Command cmd){
    // // start and stop positions read from input line
    // long start, stop;
    // Bound bound;
    // // Name of query input (e.g. AT1G01010)
    // char seqname[NAME_BUFFER_SIZE];
    // // Index of query chromosome
    // char contig_seqname[NAME_BUFFER_SIZE];
    // // query contig
    // Contig* qcon;
    // 
    // char *line = (char *) malloc(LINE_BUFFER_SIZE * sizeof(char));
    // while (fgets(line, LINE_BUFFER_SIZE, intfile) && !feof(intfile)) {
    // 
    //     // skip comments
    //     if (line[0] == '#')
    //         continue;
    // 
    //     if (!sscanf(line,
    //                 "%s %*s %*s %zu %zu %*s %*c %*s %s\n",
    //                 contig_seqname, &start, &stop, seqname)) {
    //         printf("invalid input\n");
    //         exit(EXIT_FAILURE);
    //     }
    //     qcon = syn->get_contig(0, contig_seqname);
    // 
    //     if(qcon == NULL) {
    //         fprintf(stderr, "SKIPPING ENTRY: Synteny map has no contig names '%s'\n", contig_seqname);
    //         continue;
    //     }
    // 
    //     check_in_offset(start, stop);
    //     start -= Offsets::in_start;
    //     stop -= Offsets::in_stop;
    // 
    //     bound.start = start;
    //     bound.stop = stop;
    // 
    //     switch(cmd){
    //         case C_FILTER:
    //             fprintf(stderr, "Filter function currently unavailable\n");
    //             return false;
    //         case C_COUNT:
    //             qcon->count_overlaps(bound, seqname);
    //         case C_SEARCH:
    //             qcon->find_search_intervals(bound, seqname);
    //         case C_MAP:
    //             qcon->map(bound, seqname);
    //         case C_UNSET:
    //             fprintf(stderr, "Please, give me a command\n");
    //             return false;
    //         default:
    //             fprintf(stderr, "Command '%s' not recognized\n", args.cmd.c_str());
    //             args.print_help();
    //             return false;
    //     }
    // }
    // free line;
    return true;
}
