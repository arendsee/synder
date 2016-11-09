#include "synmap.h"

Synmap::Synmap(
    std::string t_synfile,
    std::string t_tclfile,
    std::string t_qclfile,
    bool   t_swap,
    int    t_k,
    double t_r,
    char   t_trans,
    std::vector<int> t_offsets
)
    :
    synfile(t_synfile),
    tclfile(t_tclfile),
    qclfile(t_qclfile),
    swap(t_swap),
    k(t_k),
    r(t_r),
    trans(t_trans)
{
    if(t_offsets.size() != offsets.size()) {
        Rcpp::stop("Offsets must be an integer vector of 4 elements");
    }
    for(int i = 0; i < offsets.size(); i++) {
        offsets[i] = t_offsets[i];
    }
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

    std::ifstream fh(synfile);

    genome[0] = new Genome("Q");
    genome[1] = new Genome("T");

    // Contig name
    std::array< std::string, 2 > seqid;
    double score;
    char   strand;
    long   start[2], stop[2];

    Block *qblk, *tblk;

    size_t i = swap ? 1 : 0;
    size_t j = swap ? 0 : 1;

    std::string line;
    while (std::getline(fh, line)) {

        // skip comments
        if (line[0] == '#')
            continue;

        std::istringstream row(line);

        row >> seqid[i] >> start[i] >> stop[i]
            >> seqid[j] >> start[j] >> stop[j]
            >> score >> strand;

        start[0] -= offsets[0];
        start[1] -= offsets[0];
        stop[0]  -= offsets[1];
        stop[1]  -= offsets[1];

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

Rcpp::CharacterVector Synmap::filter(std::string intfile)
{

    std::ifstream fh(intfile);

    if(! fh){
        Rcpp::stop("Failed to open filter file\n");
    }

    std::string qseqid, tseqid;
    long qstart, qstop, tstart, tstop;

    Rcpp::CharacterVector out;

    std::string line;
    while (std::getline(fh, line)) {

        // skip comments
        if (line[0] == '#')
            continue;

        std::stringstream row(line);

        if (
            row >> qseqid >> qstart >> qstop
                >> tseqid >> tstart >> tstop
        ) {

            qstart -= offsets[0];
            tstart -= offsets[0];
            qstop  -= offsets[1];
            tstop  -= offsets[1];

            Feature qfeat(qseqid.c_str(), qstart, qstop);
            Feature tfeat(tseqid.c_str(), tstart, tstop);

            Contig* qcon = get_contig(0, qseqid.c_str());
            if(qcon == nullptr) {
                // Absence of a particular contig does not necessarily imply bad
                // input. So no need to throw an exception.
                Rcpp::warning("Contig '" + qseqid + "' not found in synteny map, skipping");
            } else {
                std::vector<SearchInterval> si = qcon->list_search_intervals(qfeat, r);
                for(auto &s : si) {
                    if(s.feature_overlap(&tfeat)) {
                        out.push_back(line);
                        break;
                    }
                }
            }

        } else {
            Rcpp::warning("Failed to parse line:\n\t" + line);
        }
    }

    return out;
}


std::vector<Feature> Synmap::gff2features(std::string gfffile)
{

    std::ifstream fh(gfffile);

    if(! fh){
        Rcpp::stop("Failed to open GFF file\n");
    }

    // start and stop positions read from input line
    long start, stop;
    // Name of query input (e.g. AT1G01010)
    std::string seqname;
    // Index of query chromosome
    std::string contig_seqname;
    // query contig
    Contig* qcon;

    std::vector<Feature> feats;

    std::string line;
    while (std::getline(fh, line)) {

        // skip comments
        if (line[0] == '#')
            continue;

        std::stringstream row(line);

        std::string r2, r3, r6, r7, r8;

        if (
            row >> contig_seqname >> r2 >> r3
                >> start >> stop
                >> r6 >> r7 >> r8
                >> seqname
        ){
            // check_in_offset(start, stop);
            start -= offsets[2];
            stop  -= offsets[3];
            qcon = get_contig(0, contig_seqname.c_str());
            if(qcon == nullptr) {
                char* msg;
                Rcpp::warning(
                    "SKIPPING ENTRY: Synteny map has no contig named '" +
                    std::string(contig_seqname) +
                    "'\n"
                );
                continue;
            }
            Feature feat(contig_seqname.c_str(), start, stop, seqname.c_str(), 0);
            feats.push_back(feat);
        } else {
            Rcpp::warning("Failed to parse line:\n\t" + line);
        }
    }

    return feats;
}

Rcpp::DataFrame Synmap::count(std::string intfile)
{

    std::vector<Feature> feats = gff2features(intfile);

    CountType out;

    for(auto &feat : feats) {

        Contig* qcon = get_contig(0, feat.parent_name.c_str());

        qcon->count(feat, out);
    }

    return out.as_data_frame();

}

Rcpp::DataFrame Synmap::map(std::string intfile)
{

    std::vector<Feature> feats = gff2features(intfile);

    MapType out;

    for(auto &feat : feats) {

        Contig* qcon = get_contig(0, feat.parent_name.c_str());

        qcon->map(feat, out);
    }

    return out.as_data_frame();

}

Rcpp::DataFrame Synmap::search(std::string intfile)
{

    std::vector<Feature> feats = gff2features(intfile);

    SIType out;

    for(auto &feat : feats) {

        Contig* qcon = get_contig(0, feat.parent_name.c_str());

        qcon->find_search_intervals(feat, r, out);
    }

    return out.as_data_frame();

}
