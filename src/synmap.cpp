#include "synmap.h"

Synmap::Synmap(Arguments& args)
    :
    synfile(args.synfile),
    tclfile(args.tclfile),
    qclfile(args.qclfile),
    swap(args.swap),
    k(args.k),
    trans(args.trans)
{
    load_blocks();
    validate();
}


Synmap::~Synmap()
{
    delete genome[0];
    delete genome[1];
}


void Synmap::load_blocks()
{
    if (synfile == nullptr) {
        fprintf(stderr, "NULL synteny input file (%s:%d in %s)\n", __FILE__, __LINE__, __func__);
        exit(EXIT_FAILURE);
    }

    genome[0] = new Genome("Q");
    genome[1] = new Genome("T");

    // read loop variables
    int     status = 0;
    char*   line   = nullptr;
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
                // l := -log(S) (e-values or p-values)\n"
                score = -1 * log(score);
                break;
            case 'd':
                // d := L * S (score densities)\n"
                score = score * std::min((stop[0] - start[0] + 1), (stop[1] - start[1] + 1));
                break;
            case 'p':
                // p := L * S / 100 (percent identity)\n"
                score = score * std::min((stop[0] - start[0] + 1), (stop[1] - start[1] + 1)) / 100.0;
                break;
            case 'i':
                // i := S  (default, no transformation)\n"
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
        LinkedInterval<Block>::link_homologs(qblk, tblk);

    }
    free(line);

    link_blocks();
}

void Synmap::link_blocks()
{
    genome[0]->set_contig_lengths(qclfile);
    genome[1]->set_contig_lengths(tclfile);

    genome[0]->link_block_corners();
    genome[1]->link_block_corners();

    genome[0]->set_contig_corners();
    genome[1]->set_contig_corners();

    long offset = 0;
    genome[0]->set_overlap_group(offset);
    genome[1]->set_overlap_group(offset);

    genome[0]->merge_overlaps();
    genome[0]->refresh();
    genome[1]->refresh();

    genome[0]->link_adjacent_blocks();
    genome[1]->link_adjacent_blocks();

    size_t setid = 0;
    genome[0]->link_contiguous_blocks(k, setid);
    genome[0]->transfer_contiguous_sets(genome[1]);
}


Contig* Synmap::get_contig(size_t gid, char* contig_name)
{
    if (gid == 0 || gid == 1) {
        return genome[gid]->get_contig(contig_name);
    } else {
        return nullptr;
    }
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
    genome[0]->dump_blocks();
}

void Synmap::validate()
{
    genome[0]->validate();
    genome[1]->validate();
}

void Synmap::filter(FILE* intfile)
{
    // read loop variables
    char*   line   = nullptr;
    size_t  len    = 0;
    ssize_t read;

    // Contig name
    char qseqid[NAME_BUFFER_SIZE];
    char tseqid[NAME_BUFFER_SIZE];
    char* seqid[2] = {qseqid, tseqid};
    long start[2], stop[2];

    while ((read = getline(&line, &len, intfile)) != EOF) {

        sscanf(line, "%s %ld %ld %s %ld %ld",
               seqid[0], &start[0], &stop[0],
               seqid[1], &start[1], &stop[1]);

        Feature qfeat(seqid[0], start[0], stop[0]);
        Feature tfeat(seqid[1], start[1], stop[1]);

        Contig* qcon = get_contig(0, seqid[0]);
        if(qcon == nullptr){
            // Absence of a particular contig does not necessarily imply bad
            // input. So no need to throw an exception.
            std::cerr << "WARNING: Contig '" << seqid[0] << "' not found in synteny map, skipping\n";
        } else {
            std::vector<SearchInterval> si = qcon->list_search_intervals(qfeat);

            for(auto &s : si){
                if(s.feature_overlap(&tfeat)){
                    printf("%s", line);
                    break;
                }
            }
        }
    }
    free(line);
}

bool Synmap::process_gff(FILE* intfile, Command cmd)
{
    // start and stop positions read from input line
    long start, stop;
    // Name of query input (e.g. AT1G01010)
    char seqname[NAME_BUFFER_SIZE];
    // Index of query chromosome
    char contig_seqname[NAME_BUFFER_SIZE];
    // query contig
    Contig* qcon;

    char *line = (char *) malloc(LINE_BUFFER_SIZE * sizeof(char));
    while (fgets(line, LINE_BUFFER_SIZE, intfile) && !feof(intfile)) {

        // skip comments
        if (line[0] == '#')
            continue;

        if (!sscanf(line,
                    "%s %*s %*s %ld %ld %*s %*c %*s %s\n",
                    contig_seqname, &start, &stop, seqname))
        {
            throw "invalid input in Synmap::process_gff";
        }

        // check_in_offset(start, stop);
        start -= Offsets::in_start;
        stop -= Offsets::in_stop;

        qcon = get_contig(0, contig_seqname);

        Feature feat(contig_seqname, start, stop, seqname, 0);

        if(qcon == nullptr) {
            fprintf(stderr, "SKIPPING ENTRY: Synteny map has no contig names '%s'\n", contig_seqname);
            continue;
        }

        switch(cmd){
            case C_FILTER:
                fprintf(stderr, "Filter function currently unavailable\n");
                return false;
            case C_COUNT:
                qcon->count(feat);
                break;
            case C_SEARCH:
                qcon->find_search_intervals(feat);
                break;
            case C_MAP:
                qcon->map(feat);
                break;
            case C_UNSET:
                fprintf(stderr, "Please, give me a command\n");
                return false;
            default:
                fprintf(stderr, "Command not recognized\n");
                return false;
        }
    }
    free(line);
    return true;
}
