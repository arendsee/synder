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

test_that(
  "make_tblastn_gene_map",
  {
    expect_equal(
      names(make_tblastn_gene_map(tblastn, algff)), 
      c("fprot_id", "fprot_start", "fprot_stop", "tchr", "tn_start", "tn_stop",
        "tn_evalue", "tn_ppos", "tn_frame", "tfeat_name", "tfeat_start",
        "tfeat_stop", "tfeat_src", "tfeat_type", "tfeat_strand")
    )
  }
)

# test_that(
#   "make_tblastn_si_map",
#   {
#     synres <- synder::search(tgff[1000*(1:10), ], algff)
#     expect_equal(
#       make_tblastn_si_map(tblastn, synres)
#     )
#   }
# )
