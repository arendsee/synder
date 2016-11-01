#ifndef __TYPES_H__
#define __TYPES_H__

#include <vector>
#include <Rcpp.h>

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
            Rcpp::Named("seqname") = seqname,
            Rcpp::Named("count") = count
        );
    }
};

class MapType {
private:
    std::vector<std::string> qname;
    std::vector<std::string> qchr;
    std::vector<long>        qstart;
    std::vector<long>        qstop;
    std::vector<std::string> tchr;
    std::vector<long>        tstart;
    std::vector<long>        tstop;
    std::vector<char>        strand;
    std::vector<bool>        missing;

public:
    void add_row(
        std::string t_qname,
        std::string t_qchr,
        long        t_qstart,
        long        t_qstop,
        std::string t_tchr,
        long        t_tstart,
        long        t_tstop,
        char        t_strand,
        bool        t_missing
    )
    {
        qname.push_back   ( t_qname   );
        qchr.push_back    ( t_qchr    );
        qstart.push_back  ( t_qstart  );
        qstop.push_back   ( t_qstop   );
        tchr.push_back    ( t_tchr    );
        tstart.push_back  ( t_tstart  );
        tstop.push_back   ( t_tstop   );
        strand.push_back  ( t_strand  );
        missing.push_back ( t_missing );
    }

    Rcpp::DataFrame as_data_frame() {
        return Rcpp::DataFrame::create(
            Rcpp::Named("qname")   = qname,
            Rcpp::Named("qchr")    = qchr,
            Rcpp::Named("qstart")  = qstart,
            Rcpp::Named("qstop")   = qstop,
            Rcpp::Named("tchr")    = tchr,
            Rcpp::Named("tstart")  = tstart,
            Rcpp::Named("tstop")   = tstop,
            Rcpp::Named("strand")  = strand,
            Rcpp::Named("missing") = missing
        );
    }
};

class SIType {
private:
    std::vector<std::string> name;
    std::vector<std::string> qchr;
    std::vector<long>        qstart;
    std::vector<long>        qstop;
    std::vector<std::string> tchr;
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
        std::string t_name,
        std::string t_qchr,
        long        t_qstart,
        long        t_qstop,
        std::string t_tchr,
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
        name.push_back      ( t_name      );
        qchr.push_back      ( t_qchr      );
        qstart.push_back    ( t_qstart    );
        qstop.push_back     ( t_qstop     );
        tchr.push_back      ( t_tchr      );
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
            Rcpp::Named("name")      = name,
            Rcpp::Named("qchr")      = qchr,
            Rcpp::Named("qstart")    = qstart,
            Rcpp::Named("qstop")     = qstop,
            Rcpp::Named("tchr")      = tchr,
            Rcpp::Named("tstart")    = tstart,
            Rcpp::Named("tstop")     = tstop,
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
