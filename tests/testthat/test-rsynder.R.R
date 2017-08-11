context("rsynder.R")

# ------------------------------------------------------------------------------
# Utility functions and setup
# ------------------------------------------------------------------------------

OFFSET=c(0L,0L,0L,0L,0L,0L)

df_equal <- function(obs, exp, skip=NULL){
  obs <- obs[do.call(order, as.list(obs)),]
  exp <- exp[do.call(order, as.list(exp)),]
  if(all(dim(obs) != dim(exp))){
    are_equal = FALSE
  } else {
    indices = 1:nrow(obs)
    if(!is.null(skip)){
      obs <- obs[-skip]
      exp <- exp[-skip]
    }
    are_equal = (obs == exp) %>% lapply(all) %>% unlist %>% all
  }
  are_equal
}

get_obs_exp <- function(dir, base, offsets=OFFSET, ext="", add_tcl=FALSE, add_qcl=FALSE, ...){
    syn_file <- file.path(dir, 'map.syn')
    gff_file <- file.path(dir, paste0(base, '.gff'))
    exp_file <- file.path(dir, paste0(base, ext, '-exp.txt'))
    tcl_file <- if(add_tcl) file.path(dir, 'tgen.tab') else ""
    qcl_file <- if(add_qcl) file.path(dir, 'qgen.tab') else ""

    obs <- synder::search(
      syn     = syn_file,
      gff     = gff_file,
      offsets = offsets,
      tcl     = tcl_file,
      qcl     = qcl_file,
      ...
    ) %>% as.data.frame 
    exp <- readr::read_tsv(exp_file, col_names=FALSE)
    list(obs=obs, exp=exp)
}

compare_factory <- function(dir){
  function(base, skip=c(9,10), ...){
      res <- get_obs_exp(dir, base, ...)
      df_equal(res$obs, res$exp, skip=skip)
  }
}


# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------

test_that(
  "Test IO issues",
  {
    comp <- compare_factory('io')
    expect(comp('long'), 'long GFF 9th column entries') 
  }
)

test_that(
  "Mappings beyond the edges of target scaffold",
  {
    comp <- compare_factory("unassembled")
    expect(comp('lo'),     "query is below scaffold")
    expect(comp('adj-lo'), "query is just below the scaffold")
    expect(comp('adj-hi', add_tcl=TRUE), "query is just above the scaffold")
    expect(comp('hi',     add_tcl=TRUE), "query is above the scaffold")
    expect(comp('lo', offsets=c(0L,0L,0L,0L,1L,1L), ext='-o000011'), "test with 1-base")
  }
)

test_that(
  "Test multi-chromosome cases for varying k",
  {
    #  T   =====[-------------]=====
    #        |                   |
    #  Q   ===== <->   =====   =====
    #                    |
    #  T        ===[--]=====
    #            |
    #  Q        ===
    comp <- compare_factory("interruptions/one-query-side")
    expect(comp('beside', k=0L), "query side")

    #           [----------------]
    # T             ===    ===
    #                |      |
    # Q    =====    ===    ===    =====
    #        |                      |
    # T    =====   <-->           =====
    comp <- compare_factory("interruptions/two-target-side")
    expect(comp('beside', k=0L), "target side, k=0")
    expect(comp('beside', k=1L), "target side, k=1 (should be same k=2)")

    # T    =====[--------------------]=====
    #        |                          |
    # Q    =====   ===== <--> =====   =====
    #                |          |
    # T            =====[----]=====
    comp <- compare_factory("interruptions/two-query-side")
    expect(comp('between', k=0L), "between two interior query-side intervals (k=0)")

    # T    =====[------------------------------------]=====
    #        |                                          |
    # Q    =====   =====   ===== <--> =====   =====   =====
    #                |       |          |       |
    # T            =====[----|----------|----]=====
    #                        |          |
    #                      =====[----]=====
    comp <- compare_factory("interruptions/nested")
    expect(comp('between', k=4L, ext='-k4'), "query nested two pairs of interrupting intervals (k=4)")
    expect(comp('between', k=3L, ext='-k3'), "query nested two pairs of interrupting intervals (k=3)")
  }
)

# test_that(
#   "Long 9th column in GFF",
#   {
#     # IO test, for Synder issue #2
#     synmap <- "io/map.syn"
#     # io_test "$g_synder -i $g_dir/a.gff -s $g_dir/map.syn" "long 9th GFF column"
#   }
# )

test_that(
  "One-block synteny map corner cases is correct",
  {
    comp <- compare_factory('one-block')
    expect(comp('hi'), "Query upstream of block")
    expect(comp('in'), "Query in block")
    expect(comp('lo'), "Query downstream of block")
  }
)

test_that(
  "Linear two-block synteny map is correct",
  {
    comp <- compare_factory('two-block')
    expect(comp('hi'),      "Query downstream of all blocks" )
    expect(comp('between'), "Query between two blocks"       )
    expect(comp('lo'),      "Query upstream of all blocks"   )
  }
)

test_that(
  "Linear multi-block synteny map is correct",
  {
    comp <- compare_factory('multi-block')
    expect(comp("a"), "extreme left"                              )
    expect(comp("b"), "inbetween two adjacent blocks"             )
    expect(comp("c"), "starts inbetween adjacent blocks"          )
    expect(comp("d"), "stops inbetween adjacent blocks"           )
    expect(comp("e"), "inbetween two adjacent blocks"             )
    expect(comp("f"), "starts before block 3, ends after block 3" )
    expect(comp("g"), "starts in block 2, ends after block 3"     )
    expect(comp("h"), "starts before block 2, ends after block 3" )
    expect(comp("i"), "starts in block 2, ends in block 2"        )
    expect(comp("j"), "extreme right"                             )
  }
)

test_that(
  "Simple tandem duplication",
  {
    comp <- compare_factory('simple-duplication')
    expect(comp('between'), "Query starts between the duplicated intervals")
  }
)

test_that(
  "Test when a single interval is inverted",
  {
    comp <- compare_factory('one-interval-inversion')
    expect(comp('between'), "query next to inverted interval"  )
    expect(comp('over'),    "query overlaps inverted interval" )
  }
)

test_that(
  "Test when two interval are inverted",
  {
    comp <- compare_factory("two-interval-inversion")
    expect(comp('beside'),   "query next to inverted interval"  )
    expect(comp('within'),   "query between inverted intervals" )
    expect(comp('spanning'), "query spans inverted intervals"   )
  }
)

test_that(
  "Test tandem transposition",
  {
    comp <- compare_factory("tandem-transposition")
    expect(comp('beside'), "query beside the transposed pair")
    expect(comp('within'), "query between the transposed pair")
  }
)

test_that(
  "Test target side internal overlaps",
  {
    comp <- compare_factory("irregular-overlaps")
    # This can fail if you are 1) not sorting the by_stop vector in Contig by
    # Block stop positions, or 2) are snapping the search interval left
    # boundary to a Block that is nearest by start, but not be stop.
    expect(comp('left'), "left side")
    expect(comp('right'), "right side")
  }
)

test_that(
  "Test overlap edge cases",
  {
    comp <- compare_factory("off-by-one")
    expect(comp('a'), "overlap of 1")
  }
)

test_that(
  "Extreme value resulting from an inversion",
  {
    comp <- compare_factory("inverted-extremes")
    expect(comp('extreme'), "between the query intervals, extreme SI")
  }
)

test_that(
  "Deletion tests (adjacent bounds in target)",
  {
    comp <- compare_factory("deletion")
    expect(comp('between'), "query is inbetween")
  }
)

test_that(
  "Test multi-chromosome cases when k=0",
  {
    comp_mc <- compare_factory("interruptions/multi-chromosome")
    #  T   =====[---->
    #        |
    #  Q   =====   <->   =====
    #                      |
    #  T           <----]=====
    expect(comp_mc('between'), "interuption between query intervals")

    # T             ===    ===
    #                |      |
    # Q    =====[--]===    ===[--]=====
    #        |                      |
    # T    =====   <-->           =====
    comp_1q <- compare_factory("interruptions/one-query-side")
    expect(comp_1q('beside'), "target side")

    #  T   =====[-------------]=====
    #        |                   |
    #  Q   ===== <->   =====   =====
    #                    |
    #  T        ===[--]=====
    #            |
    #  Q        ===
    expect(comp_1q('beside'), "query side")


    # T    =====                      =====
    #        |                          |
    # Q    =====   ===== <--> =====   =====
    #                |          |
    # T            =====[----]=====
    comp_2q <- compare_factory("interruptions/two-query-side")
    expect(comp_2q('between', k=0L), "between two interior query-side intervals (k=0)")

  }
)


test_that(
  "Confirm two-scaffold systems are unaffected by k",
  {
    comp <- compare_factory("tandem-transposition")
    expect(comp('beside', k=4L), "query beside the transposed pair")
    expect(comp('within', k=4L), "query between the transposed pair")

    comp <- compare_factory("simple-duplication")
    expect(comp('between', k=4L), "query starts between the duplicated intervals")
  }
)

test_that(
  "Assert overlapping regions are correctly merged",
  {
    # The `overlap` set of synteny maps are ones that require synder merge doubly
    # overlapping intervals. There are a variety of tricky corner cases. Apart from
    # the first few, all the test below were added to test for a specific bug.
    # Basically, don't mess with these.
    o1  <- synder::dump('overlap/map-1.syn',  trans='p', offsets=OFFSET) %>% as.data.frame
    o2  <- synder::dump('overlap/map-2.syn',  trans='p', offsets=OFFSET) %>% as.data.frame
    o3  <- synder::dump('overlap/map-3.syn',  trans='p', offsets=OFFSET) %>% as.data.frame
    o4  <- synder::dump('overlap/map-4.syn',  trans='p', offsets=OFFSET) %>% as.data.frame
    o5  <- synder::dump('overlap/map-5.syn',  trans='p', offsets=OFFSET) %>% as.data.frame
    o6  <- synder::dump('overlap/map-6.syn',  trans='p', offsets=OFFSET) %>% as.data.frame
    o7  <- synder::dump('overlap/map-7.syn',  trans='p', offsets=OFFSET) %>% as.data.frame
    o8  <- synder::dump('overlap/map-8.syn',  trans='p', offsets=OFFSET) %>% as.data.frame
    o9  <- synder::dump('overlap/map-9.syn',  trans='p', offsets=OFFSET) %>% as.data.frame
    o10 <- synder::dump('overlap/map-10.syn', trans='p', offsets=OFFSET) %>% as.data.frame
    o11 <- synder::dump('overlap/map-11.syn', trans='p', offsets=OFFSET) %>% as.data.frame

    # single nesting
    expect_equal(o1[[2]], c(100,400,600))
    expect_equal(o1[[5]], c(1100,1500,1700))
    expect_equal(o1[[7]], c(101,101,101))

    # triple identical
    expect_equal(o2[[2]], c(100,400,600))
    expect_equal(o2[[5]], c(1100,1500,1700))
    expect_equal(o2[[7]], c(101,101,101))

    # left
    expect_equal(o3[[2]], c(100,400,600))
    expect_equal(o3[[5]], c(1100,1500,1700))
    expect_equal(o3[[7]], c(101,101,101))

    # left-right
    expect_equal(o4[[2]], c(100,400,600))
    expect_equal(o4[[5]], c(1100,1500,1700))
    expect_equal(o4[[7]], c(101,101,101))

    # double nest
    expect_equal(o5[[2]], c(100,400,600))
    expect_equal(o5[[5]], c(1100,1500,1700))
    expect_equal(round(o5[[7]]), c(101,141,101))

    # double nest left
    expect_equal(o6[[2]], c(100,400,600))
    expect_equal(o6[[5]], c(1100,1500,1700))
    expect_equal(round(o6[[7]]), c(101,141,101))

    # Q-inside T-right
    expect_equal(o7[[2]], c(100,400))
    expect_equal(o7[[5]], c(1100,1500))
    expect_equal(round(o7[[7]]), c(101,91))

    # Tangles
    expect_equal(o8[[2]], c(100,1100))
    expect_equal(o8[[5]], c(1100,1500))
    expect_equal(round(o8[[7]]), c(101,700))

    # double overlap on different target contigs
    expect_equal(o9[[2]], c(100,400,600,600))
    expect_equal(o9[[5]], c(1100,1500,600,600))
    expect_equal(o9[[7]], c(101,101,101,101))

    # transitive group ids
    #  T                ===C===
    #      ===A===   ==B==  |
    #         __\ _____/    |
    #        /   \          |
    #  Q   ==a==  \  ====c====
    #        ======b=====
    # overlapping groups: (abc), (A), (BC)
    # ((BC), (ac)) should NOT be merged
    expect_equal(o10[[2]], c(100,100,300))
    expect_equal(o10[[5]], c(300,100,400))
    expect_equal(round(o10[[7]]), c(101,101,201))

    #  T   =======A=======
    #          ===B===  |
    #          /        |
    #  Q      /  ====a===
    #       ===b====
    # This caused problems when A was merged into B, this left ManyBlocks::front
    # pointing to a broken Block
    expect_equal(o11[[2]], c(100))
    expect_equal(o11[[5]], c(100))
    expect_equal(round(o11[[7]]), c(251))
  }
)
