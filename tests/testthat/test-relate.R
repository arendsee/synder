context("relate")

algff_file <- system.file("arabidopsis", "al.gff", package="synder")
algff <- synder::as_gff(algff_file)

atgff_file <- system.file("arabidopsis", "at.gff", package="synder")
atgff <- synder::as_gff(atgff_file)

syn_file <- system.file("arabidopsis", "at-vs-al.syn", package="synder")
syn <- synder::as_synmap(syn_file)

blastp_file <- system.file("arabidopsis", "nfyc_at-vs-al_blastp.tab", package="synder")
blastp <- read.table(blastp_file, header=TRUE)

tblastn_file <- system.file("arabidopsis", "nfyc_at-vs-al_tblastn.tab", package="synder")
tblastn <- read.table(tblastn_file, header=TRUE)

srcres <- synder::search(syn_file, atgff)

has_names <- function(d, xs){
  all(xs %in% names(d))
}

has_prefix <- function(xs, prefix){
  pattern <- paste0("^", prefix)
  all(grepl(pattern, xs, perl=TRUE))
}

extractNameFromAttr <- function(x){
  sub(".*Name=([^;]+).*", "\\1", x)
}

test_that("tidy gff",
  {
    expect_equal(names(.tidy_gff(atgff)), c("seqid", "chr", "source", "type", "start", "stop", "strand"))
    expect_equal(names(.tidy_gff(atgff, prefix="f")), c("fseqid", "fchr", "fsource", "ftype", "fstart", "fstop", "fstrand"))
    expect_equal(unname(sapply(.tidy_gff(atgff), class)), c("character", "character", "character", "character", "integer", "integer", "character"))
  }
)

test_that("tidy search result",
  {
    rcol <- c("fseqid", "fchr", "si_fstart", "si_fstop", "tchr", "si_tstart",
              "si_tstop", "orientation", "si_score", "cset", "l_flag",
              "r_flag", "inbetween")
    expect_true(has_names(.tidy_searchResult(srcres), rcol))
  }
)

test_that("load_blastp_file",
  {
    x <- load_blastp_file(blastp_file)
    # required columns
    rcols <- c("fseqid", "tseqid", "fprotein", "tprotein", "bp_evalue", "bp_id")
    # there may be more columns than this
    expect_true(has_names(x, rcols))
    # all non-id columns have the bp_ prefix
    expect_true(has_prefix(names(x)[5:ncol(x)], "bp_"))
  }
)

test_that("load_tblastn_file",
  {
    x <- load_tblastn_file(tblastn_file)
    # required columns
    rcols <- c("fseqid", "fprotein", "tchr", "tn_evalue", "tn_sframe", "tn_id")
    # there may be more columns than this
    expect_true(has_names(x, rcols))
    # all non-id columns have the bp_ prefix
    expect_true(has_prefix(names(x)[4:ncol(x)], "tn_"))
    # tn_sframe is correctly typed
    expect_true(all(x$tn_sframe %in% c(-3,-2,-1,1,2,3)))
  }
)

test_that(
  "make_tblastn_gene_map",
  {
    x <- make_tblastn_gene_map(tblastn_file, algff)
    rcol <- c("fseqid", "fprotein", "prot_fstart", "prot_fstop", "tchr",
              "tn_tstart", "tn_tstop", "tn_tframe", "tn_evalue", "tn_ppos",
              "feat_tname", "feat_tstart", "feat_tstop", "feat_tsrc",
              "feat_ttype", "feat_tstrand")
    expect_true(has_names(x, rcol))
  }
)

test_that(
  "featureMap",
  {
    x <- featureMap(srcres, tgff=algff, fgff=atgff)
    rcol <- c("fseqid", "ftype", "fsource", "fchr", "fstart", "fstop",
              "fstrand", "tseqid", "ttype", "tsource", "tchr", "tstart",
              "tstop", "tstrand", "tstrand", "si_fstart", "si_fstop",
              "si_tstart", "si_tstop", "orientation", "si_score", "cset",
              "l_flag", "r_flag", "inbetween", "strands_agree")
    expect_true(has_names(x, rcol))
  }
)

test_that(
  "make_tblastn_si_map",
  {
    x <- make_tblastn_si_map(tblastn_file, srcres, fmap=NULL)
    rcol <- c("fseqid", "fprotein", "prot_fstart", "prot_fstop", "tchr",
              "tn_tstart", "tn_tstop", "tn_tframe", "tn_evalue", "tn_ppos",
              "si_fseqid", "cset", "fchr", "si_fstart", "si_fstop",
              "orientation", "si_score", "l_flag", "r_flag", "inbetween")
    expect_true(has_names(x, rcol))
  }
)
