#include "synmap.h"

Synmap::Synmap(
    FILE*  synfile,
    FILE*  tclfile,
    FILE*  qclfile,
    bool   swap,
    int    k,
    double r,
    char   trans
)
    :
    synfile(synfile),
    tclfile(tclfile),
    qclfile(qclfile),
    swap(swap),
    k(k),
    r(r),
    trans(trans)
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

        start[0] -= Offsets::syn_start;
        start[1] -= Offsets::syn_start;
        stop[0]  -= Offsets::syn_stop;
        stop[1]  -= Offsets::syn_stop;

        if(status != 8) {
            Rcpp::stop("Failed to read input line:\n" + std::string(line));
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
                Rcpp::stop("Unexpected value of transform (trans argument)");
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

Rcpp::DataFrame Synmap::as_data_frame()
{
    return genome[0]->as_data_frame();
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


Contig* Synmap::get_contig(size_t gid, const char* contig_name)
{
    if (gid == 0 || gid == 1) {
        return genome[gid]->get_contig(contig_name);
    } else {
        return nullptr;
    }
}

void Synmap::validate()
{
    genome[0]->validate();
    genome[1]->validate();
}

Rcpp::CharacterVector Synmap::filter(FILE* intfile)
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

    Rcpp::CharacterVector out;

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
            Rcpp::warning("Contig '" + std::string(seqid[0]) + "' not found in synteny map, skipping");
        } else {
            std::vector<SearchInterval> si = qcon->list_search_intervals(qfeat, r);

            for(auto &s : si){
                if(s.feature_overlap(&tfeat)){
                    out.push_back(line);
                    break;
                }
            }
        }
    }
    free(line);

    return out;
}


std::vector<Feature> Synmap::gff2features(FILE* fh){
    // start and stop positions read from input line
    long start, stop;
    // Name of query input (e.g. AT1G01010)
    char seqname[NAME_BUFFER_SIZE];
    // Index of query chromosome
    char contig_seqname[NAME_BUFFER_SIZE];
    // query contig
    Contig* qcon;

    std::vector<Feature> feats;

    char *line = (char *) malloc(LINE_BUFFER_SIZE * sizeof(char));
    while (fgets(line, LINE_BUFFER_SIZE, fh) && !feof(fh)) {

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

        if(qcon == nullptr) {
            std::cerr
                << "SKIPPING ENTRY: Synteny map has no contig named '"
                << contig_seqname
                << "'\n";
            continue;
        }

        Feature feat(contig_seqname, start, stop, seqname, 0);

        feats.push_back(feat);

    }
    free(line);

    return feats;
}

Rcpp::DataFrame Synmap::count(FILE* intfile)
{

    std::vector<Feature> feats = gff2features(intfile);

    CountType out;

    for(auto &feat : feats) {

        Contig* qcon = get_contig(0, feat.parent_name.c_str());

        qcon->count(feat, out);
    }

    return out.as_data_frame();

}

Rcpp::DataFrame Synmap::map(FILE* intfile)
{

    std::vector<Feature> feats = gff2features(intfile);

    MapType out;

    for(auto &feat : feats) {

        Contig* qcon = get_contig(0, feat.parent_name.c_str());

        qcon->map(feat, out);
    }

    return out.as_data_frame();

}

Rcpp::DataFrame Synmap::search(FILE* intfile)
{

    std::vector<Feature> feats = gff2features(intfile);

    SIType out;

    for(auto &feat : feats) {

        Contig* qcon = get_contig(0, feat.parent_name.c_str());

        qcon->find_search_intervals(feat, r, out);
    }

    return out.as_data_frame();

}
