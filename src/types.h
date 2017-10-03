#ifndef __TYPES_H__
#define __TYPES_H__

#include <vector>
#include <Rcpp.h>

// Convert to 1-based position vector
std::vector<long> to_one_base(std::vector<long> x);

class DumpType {
private:
    std::vector<std::string> qcon;
    std::vector<long>        qstart;
    std::vector<long>        qstop;
    std::vector<std::string> tcon;
    std::vector<long>        tstart;
    std::vector<long>        tstop;
    std::vector<double>      score;
    std::vector<char>        strand;
    std::vector<size_t>      cset;

public:
    void add_row(
        std::string t_qcon,
        long        t_qstart,
        long        t_qstop,
        std::string t_tcon,
        long        t_tstart,
        long        t_tstop,
        double      t_score,
        char        t_strand,
        size_t      t_cset
    )
    {
        qcon.push_back   ( t_qcon   );
        qstart.push_back ( t_qstart );
        qstop.push_back  ( t_qstop  );
        tcon.push_back   ( t_tcon   );
        tstart.push_back ( t_tstart );
        tstop.push_back  ( t_tstop  );
        score.push_back  ( t_score  );
        strand.push_back ( t_strand );
        cset.push_back   ( t_cset   );
    }

    Rcpp::DataFrame as_data_frame() {
        return Rcpp::DataFrame::create(
            Rcpp::Named("qseqid") = qcon,
            Rcpp::Named("qstart") = to_one_base(qstart),
            Rcpp::Named("qstop")  = to_one_base(qstop),
            Rcpp::Named("tseqid") = tcon,
            Rcpp::Named("tstart") = to_one_base(tstart),
            Rcpp::Named("tstop")  = to_one_base(tstop),
            Rcpp::Named("score")  = score,
            Rcpp::Named("strand") = strand,
            Rcpp::Named("cset")   = cset
        );
    }
};

class CountType {
private:
    std::vector<std::string> seqname;
    std::vector<int> count;

public:
    void add_row(std::string s, int c) {
        seqname.push_back(s);
        count.push_back(c);
    }

    Rcpp::DataFrame as_data_frame() {
        return Rcpp::DataFrame::create(
            // NOTE: what I call the seqname internally in C synder, comes
            // from the 9th GFF column (currently), so technically shouldn't be
            // called the 'seqname'.
            Rcpp::Named("attr")  = seqname,
            Rcpp::Named("count") = count
        );
    }
};

class MapType {
private:
    std::vector<std::string> seqname;
    std::vector<std::string> qcon;
    std::vector<long>        qstart;
    std::vector<long>        qstop;
    std::vector<std::string> tcon;
    std::vector<long>        tstart;
    std::vector<long>        tstop;
    std::vector<char>        strand;
    std::vector<bool>        missing;

public:
    void add_row(
        std::string t_seqname,
        std::string t_qcon,
        long        t_qstart,
        long        t_qstop,
        std::string t_tcon,
        long        t_tstart,
        long        t_tstop,
        char        t_strand,
        bool        t_missing
    )
    {
        seqname.push_back ( t_seqname );
        qcon.push_back    ( t_qcon    );
        qstart.push_back  ( t_qstart  );
        qstop.push_back   ( t_qstop   );
        tcon.push_back    ( t_tcon    );
        tstart.push_back  ( t_tstart  );
        tstop.push_back   ( t_tstop   );
        strand.push_back  ( t_strand  );
        missing.push_back ( t_missing );
    }

    Rcpp::DataFrame as_data_frame() {
        return Rcpp::DataFrame::create(
            Rcpp::Named("attr")    = seqname,
            Rcpp::Named("qseqid")  = qcon,
            Rcpp::Named("qstart")  = to_one_base(qstart),
            Rcpp::Named("qstop")   = to_one_base(qstop),
            Rcpp::Named("tseqid")  = tcon,
            Rcpp::Named("tstart")  = to_one_base(tstart),
            Rcpp::Named("tstop")   = to_one_base(tstop),
            Rcpp::Named("strand")  = strand,
            Rcpp::Named("missing") = missing
        );
    }
};

class SIType {
private:
    std::vector<std::string> seqname;
    std::vector<std::string> qcon;
    std::vector<long>        qstart;
    std::vector<long>        qstop;
    std::vector<std::string> tcon;
    std::vector<long>        tstart;
    std::vector<long>        tstop;
    std::vector<char>        strand;
    std::vector<double>      score;
    std::vector<size_t>      cset;
    std::vector<int>         l_flag;
    std::vector<int>         r_flag;
    std::vector<bool>        inbetween;

public:
    void add_row(
        std::string t_seqname,
        std::string t_qcon,
        long        t_qstart,
        long        t_qstop,
        std::string t_tcon,
        long        t_tstart,
        long        t_tstop,
        char        t_strand,
        double      t_score,
        size_t      t_cset,
        int         t_l_flag,
        int         t_r_flag,
        bool        t_inbetween
    )
    {
        seqname.push_back   ( t_seqname   );
        qcon.push_back      ( t_qcon      );
        qstart.push_back    ( t_qstart    );
        qstop.push_back     ( t_qstop     );
        tcon.push_back      ( t_tcon      );
        tstart.push_back    ( t_tstart    );
        tstop.push_back     ( t_tstop     );
        strand.push_back    ( t_strand    );
        score.push_back     ( t_score     );
        cset.push_back      ( t_cset      );
        l_flag.push_back    ( t_l_flag    );
        r_flag.push_back    ( t_r_flag    );
        inbetween.push_back ( t_inbetween );
    }

    Rcpp::DataFrame as_data_frame() {
        return Rcpp::DataFrame::create(
            Rcpp::Named("attr")      = seqname,
            Rcpp::Named("qseqid")    = qcon,
            Rcpp::Named("qstart")    = to_one_base(qstart),
            Rcpp::Named("qstop")     = to_one_base(qstop),
            Rcpp::Named("tseqid")    = tcon,
            Rcpp::Named("tstart")    = to_one_base(tstart),
            Rcpp::Named("tstop")     = to_one_base(tstop),
            Rcpp::Named("strand")    = strand,
            Rcpp::Named("score")     = score,
            Rcpp::Named("cset")      = cset,
            Rcpp::Named("l_flag")    = l_flag,
            Rcpp::Named("r_flag")    = r_flag,
            Rcpp::Named("inbetween") = inbetween
        );
    }
};

#endif
