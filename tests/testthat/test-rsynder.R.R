context("rsynder.R")

OFFSET=c(0,0,0,0,0,0)

df_equal <- function(obs, exp, skip=NULL){
  if(ncol(obs) != ncol(exp)){
    FALSE
  } else {
    indices = 1:nrow(obs)
    if(!is.null(skip)){
      obs <- obs[-skip]
      exp <- exp[-skip]
    }
    (obs == exp) %>% lapply(all) %>% unlist %>% all
  }
}

test_that(
  "One-block synteny map corner cases is correct",
  {
    synmap <- 'one-block/map.syn'

    hi_obs <- synder::search(synmap, 'one-block/hi.gff', offsets=OFFSET)
    in_obs <- synder::search(synmap, 'one-block/in.gff', offsets=OFFSET)
    lo_obs <- synder::search(synmap, 'one-block/lo.gff', offsets=OFFSET)

    hi_exp <- readr::read_tsv('one-block/hi-exp.txt', col_names=FALSE)
    in_exp <- readr::read_tsv('one-block/in-exp.txt', col_names=FALSE)
    lo_exp <- readr::read_tsv('one-block/lo-exp.txt', col_names=FALSE)

    expect(df_equal(hi_obs, hi_exp, skip=c(9,10)), "Query upstream of block")
    expect(df_equal(in_obs, in_exp, skip=c(9,10)), "Query in block")
    expect(df_equal(lo_obs, lo_exp, skip=c(9,10)), "Query downstream of block")
  }
)

test_that(
  "Two-block synteny map is correct",
  {
    synmap <- 'two-block/map.syn'

    hi_obs <- synder::search(synmap, 'two-block/hi.gff', offsets=OFFSET)
    hi_exp <- readr::read_tsv('two-block/hi-exp.txt', col_names=FALSE)
    expect(df_equal(hi_obs, hi_exp, skip=c(9,10)), "Query downstream of all blocks")

    between_obs <- synder::search(synmap, 'two-block/between.gff', offsets=OFFSET)
    between_exp <- readr::read_tsv('two-block/between-exp.txt', col_names=FALSE)
    expect(df_equal(between_obs, between_exp, skip=c(9,10)), "Query between two blocks")

    lo_obs <- synder::search(synmap, 'two-block/lo.gff', offsets=OFFSET)
    lo_exp <- readr::read_tsv('two-block/lo-exp.txt', col_names=FALSE)
    expect(df_equal(lo_obs, lo_exp, skip=c(9,10)), "Query upstream of all blocks")
  }
)

# The `overlap` set of synteny maps are ones that require synder merge doubly
# overlapping intervals. There are a variety of tricky corner cases. Apart from
# the first few, all the test below were added to test for a specific bug.
# Basically, don't mess with these.
o1  <- synder::dump('overlap/map-1.syn',  trans='p', offsets=OFFSET)
o2  <- synder::dump('overlap/map-2.syn',  trans='p', offsets=OFFSET)
o3  <- synder::dump('overlap/map-3.syn',  trans='p', offsets=OFFSET)
o4  <- synder::dump('overlap/map-4.syn',  trans='p', offsets=OFFSET)
o5  <- synder::dump('overlap/map-5.syn',  trans='p', offsets=OFFSET)
o6  <- synder::dump('overlap/map-6.syn',  trans='p', offsets=OFFSET)
o7  <- synder::dump('overlap/map-7.syn',  trans='p', offsets=OFFSET)
o8  <- synder::dump('overlap/map-8.syn',  trans='p', offsets=OFFSET)
o9  <- synder::dump('overlap/map-9.syn',  trans='p', offsets=OFFSET)
o10 <- synder::dump('overlap/map-10.syn', trans='p', offsets=OFFSET)
o11 <- synder::dump('overlap/map-11.syn', trans='p', offsets=OFFSET)

test_that(
  "Assert overlapping regions are correctly merged",
  {
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


# - test not yet merged
# # ---------------------------------------------------------------------
# g_dir="$PWD/test/test-data/multi-block"
# announce "\nTesting with 5 adjacent blocks on the same strand"
# runtest a "extreme left"
# runtest b "inbetween two adjacent blocks"
# runtest c "starts inbetween adjacent blocks"
# runtest d "stops inbetween adjacent blocks"
# runtest e "inbetween two adjacent blocks"
# runtest f "starts before block 3, ends after block 3"
# runtest g "starts in block 2, ends after block 3"
# runtest h "starts before block 2, ends after block 3"
# runtest i "starts in block 2, ends in block 2"
# runtest j "extreme right"
#
# # ---------------------------------------------------------------------
# g_dir="$PWD/test/test-data/simple-duplication"
# announce "\nTest simple tandem duplication"
# runtest between "query starts between the duplicated intervals"
#
# # ---------------------------------------------------------------------
# g_dir="$PWD/test/test-data/one-interval-inversion"
# announce "\nTest when a single interval is inverted"
# runtest between "query next to inverted interval"
# runtest over    "query overlaps inverted interval"
#
# # ---------------------------------------------------------------------
# g_dir="$PWD/test/test-data/two-interval-inversion"
# announce "\nTest when two interval are inverted"
# runtest beside   "query next to inverted interval"
# runtest within   "query between inverted intervals"
# runtest spanning "query spans inverted intervals"
#
# # ---------------------------------------------------------------------
# g_dir="$PWD/test/test-data/tandem-transposition"
# announce "\nTest tandem transposition"
# runtest beside "query beside the transposed pair"
# runtest within "query between the transposed pair"
#
# # ---------------------------------------------------------------------
# g_dir="$PWD/test/test-data/irregular-overlaps"
# announce "\nTest target side internal overlaps"
# runtest left "left side" "You are either 1) not sorting the by_stop vector
# in Contig by Block stop positions, or 2) are snapping the search interval left
# boundary to a Block that is nearest by start, but not be stop."
# runtest right "right side"
#
# # ---------------------------------------------------------------------
# g_dir="$PWD/test/test-data/off-by-one"
# announce "\nTest overlap edge cases"
# runtest a "overlap of 1"
#
# # ---------------------------------------------------------------------
# g_dir="$PWD/test/test-data/inverted-extremes"
# announce "\nExtreme value resulting from an inversion"
# runtest extreme "between the query intervals, extreme SI"
#
# # ---------------------------------------------------------------------
# g_dir="$PWD/test/test-data/deletion"
# announce "\nDeletion tests (adjacent bounds in target)"
# runtest between "query is inbetween"
#
# # ---------------------------------------------------------------------
# g_dir="$PWD/test/test-data/unassembled"
# announce "\nMappings beyond the edges of target scaffold"
# runtest lo "query is below scaffold"
# runtest adj-lo "query is just below the scaffold"
# runtest adj-hi "query is just above the scaffold"
# runtest hi "query is above the scaffold"
# runtest lo "test with 1-base" 0 1
#
# # ---------------------------------------------------------------------
# announce "\nTest multi-chromosome cases when k=0"
# g_arg=" -k 0 "
# #  T   =====[---->
# #        |
# #  Q   =====   <->   =====
# #                      |
# #  T           <----]=====
# g_dir="$PWD/test/test-data/interruptions/multi-chromosome"
# runtest between "interuption between query intervals"
# #  T   =====[-------------]=====
# #        |                   |
# #  Q   ===== <->   =====   =====
# #                    |
# #  T        ===[--]=====
# #            |
# #  Q        ===
# g_dir="$PWD/test/test-data/interruptions/one-query-side"
# runtest beside "query side"
# # T             ===    ===
# #                |      |
# # Q    =====[--]===    ===[--]=====
# #        |                      |
# # T    =====   <-->           =====
# g_dir="$PWD/test/test-data/interruptions/two-target-side"
# runtest beside "target side"
# g_arg=" -k 1 "
# runtest beside "target side, k=1 (should be the same)"
# # T    =====                      =====
# #        |                          |
# # Q    =====   ===== <--> =====   =====
# #                |          |
# # T            =====[----]=====
# g_dir="$PWD/test/test-data/interruptions/two-query-side"
# ark=" -k 0 "
# runtest between "between two interior query-side intervals (k=0)"
#
# # ---------------------------------------------------------------------
# set_defaults
# g_arg=' -k 4 '
# announce "\nConfirm two-scaffold systems are unaffected by k"
# g_dir="$PWD/test/test-data/tandem-transposition"
# runtest beside "query beside the transposed pair"
# runtest within "query between the transposed pair"
# g_dir="$PWD/test/test-data/simple-duplication"
# runtest between "query starts between the duplicated intervals"
#
# # ---------------------------------------------------------------------
# set_defaults
# announce "\nTest multi-chromosome cases when k=2"
# g_arg=" -k 2 "
# g_exp_ext='exp-k2'
# #  T   =====[-------------]=====
# #        |                   |
# #  Q   ===== <->   =====   =====
# #                    |
# #  T        ===[--]=====
# #            |
# #  Q        ===
# g_dir="$PWD/test/test-data/interruptions/one-query-side"
# runtest beside "query side"
# #           [----------------]
# # T             ===    ===
# #                |      |
# # Q    =====    ===    ===    =====
# #        |                      |
# # T    =====   <-->           =====
# g_dir="$PWD/test/test-data/interruptions/two-target-side"
# runtest beside "target side"
# # T    =====[--------------------]=====
# #        |                          |
# # Q    =====   ===== <--> =====   =====
# #                |          |
# # T            =====[----]=====
# g_dir="$PWD/test/test-data/interruptions/two-query-side"
# runtest between "between two interior query-side intervals (k=2)"
# # T    =====[------------------------------------]=====
# #        |                                          |
# # Q    =====   =====   ===== <--> =====   =====   =====
# #                |       |          |       |
# # T            =====[----|----------|----]=====
# #                        |          |
# #                      =====[----]=====
# g_dir="$PWD/test/test-data/interruptions/nested"
# g_arg=" -k 4 "
# g_exp_ext="exp-k4"
# runtest between "query nested two pairs of interrupting intervals (k=4)"
# g_arg=" -k 3 "
# g_exp_ext="exp-k3"
# runtest between "query nested two pairs of interrupting intervals (k=3)"
#
# # ---------------------------------------------------------------------
# set_defaults
# g_dir="$PWD/test/test-data/synmap-overlaps"
# announce "\nsyntenic overlaps"
# runtest simple "Between the weird"
#
# # IO test, for Synder issue #2
# g_dir="$PWD/test/test-data/io"
# # io_test "$g_synder -i $g_dir/a.gff -s $g_dir/map.syn" "long 9th GFF column"
