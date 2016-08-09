#!/usr/bin/env bash
set -u

synder=$PWD/synder

total_passed=0
total_failed=0

announce(){
    [[ -t 1 ]] && o="\e[1;33m$1\e[0m" || o=$1
    echo -e $o
}

warn(){
    [[ -t 1 ]] && o="\e[1;31m$1\e[0m" || o=$1
    echo -en $o
}

emphasize(){
   emphasize_n "$1" 
   echo
}

emphasize_n(){
    [[ -t 1 ]] && o="\e[1;39m$1\e[0m" || o=$1
    echo -ne $o
}

# A function to select which parts of the output should be compared
filter () {
    cut -f1-7
}

runtest(){
    dif=$1
    base=$2
    msg=$3
    errmsg=${4:-0}
    tmp=/tmp/synder-$RANDOM$RANDOM
    mkdir $tmp
    echo -n "Testing $msg ... "
    $synder -d $dir/map.syn a b $tmp/db 
    $synder -i $dir/$base.gff -s $tmp/db/a_b.txt -c search > $tmp/a

    # Since flags are currently in flux, test only the first 7 columns
    diff <(cat $tmp/a | filter) \
         <(cat $dir/${base}-exp.txt | filter) > /dev/null 

    if [[ $? == 0 ]]
    then
        total_passed=$(( $total_passed + 1 ))
        echo "OK"
    else
        echo
        warn "FAIL"
        echo " (in `basename $dir`/)"
        [[ $errmsg == 0 ]] || (echo -e $errmsg | fmt)
        total_failed=$(( $total_failed + 1 ))
        echo "======================================="
        emphasize_n "expected output"; echo ": (${base}-exp.txt)"
        cat $dir/${base}-exp.txt | filter | column -t
        emphasize "observed output:"
        cat $tmp/a | filter | column -t
        emphasize_n "query gff"; echo ": (${base}.gff)"
        column -t $dir/$base.gff
        emphasize_n "synteny map"; echo ": (map.syn)"
        column -t $dir/map.syn
        echo "See zzz_g and zzz_db"
        echo -e "---------------------------------------\n"

        ln -sf $dir/$base.gff zzz_g 
        $synder -d $dir/map.syn a b zzz_db

    fi

    rm -rf $tmp
}

# #---------------------------------------------------------------------
dir="$PWD/test/test-data/one-block"
announce "\nTesting with synteny map length == 1"
runtest $dir hi     "query after of block"
runtest $dir within "query within block"
runtest $dir lo     "query before of block"

# #---------------------------------------------------------------------
dir="$PWD/test/test-data/two-block"
announce "\nTesting with synteny map length == 2"
runtest $dir hi      "query downstream of all blocks"
runtest $dir between "query between two blocks"
runtest $dir lo      "query upstream of all blocks"

# #---------------------------------------------------------------------
dir="$PWD/test/test-data/multi-block"
announce "\nTesting with 5 adjacent blocks on the same strand"
runtest $dir a "Extreme left"
runtest $dir b "Inbetween two adjacent blocks"
runtest $dir c "Starts inbetween adjacent blocks"
runtest $dir d "Stops inbetween adjacent blocks"
runtest $dir e "Inbetween two adjacent blocks"
runtest $dir f "Starts before block 3, ends after block 3"
runtest $dir g "Starts in block 2, ends after block 3"
runtest $dir h "Starts before block 2, ends after block 3"
runtest $dir i "Starts in block 2, ends in block 2"
runtest $dir j "Extreme right"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/simple-duplication"
announce "\nTest simple tandem duplication"
runtest $dir between "Query starts between the duplicated intervals"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/one-interval-inversion"
announce "\nTest when a single interval is inverted"
runtest $dir between "Query next to inverted interval"
runtest $dir over    "Query overlaps inverted interval"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/two-interval-inversion"
announce "\nTest when two interval are inverted"
runtest $dir beside   "Query next to inverted interval"
runtest $dir within   "Query between inverted intervals"
# runtest $dir spanning "Query spans inverted intervals"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/tiny-indel-query-side"
announce "\nTest when a small interval interupts on one side"
runtest $dir beside "Query side"
dir="$PWD/test/test-data/tiny-indel-target-side"
runtest $dir beside "Target side"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/tandem-transposition"
announce "\nTest tandem transposition"
runtest $dir beside "Query beside the transposed pair"
# runtest $dir within "Query between the transposed pair"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/irregular-overlaps"
announce "\nTest target side internal overlaps"
runtest $dir left "Left side" "You are either 1) not sorting the by_stop vector
in Contig by Block stop positions, or 2) are snapping the search interval left
boundary to a Block that is nearest by start, but not be stop."
# runtest $dir right "Right side"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/multi-chromosome"
announce "\nTest two intervals on same query chr but different target chr"
# runtest $dir between "Between the query intervals"


warn "And the hard ones" ; echo
    dir="$PWD/test/test-data/two-interval-inversion"
    announce "\n1: Test when two interval are inverted"
    runtest $dir spanning "Query spans inverted intervals"
    dir="$PWD/test/test-data/tandem-transposition"
    announce "\n2: Test tandem transposition"
    runtest $dir within "Query between the transposed pair"
    dir="$PWD/test/test-data/irregular-overlaps"
    announce "\n3: Test target side internal overlaps"
    runtest $dir right "Right side"

#---------------------------------------------------------------------
echo

#---------------------------------------------------------------------
# valgrind

valgrind_checked=0
valgrind_exit_status=1
if hash valgrind 2> /dev/null; then
    valgrind_checked=1
    dir="$PWD/test/test-data/multi-block"
    tmp=/tmp/synder-$RANDOM
    mkdir $tmp
    $synder -d $dir/map.syn a b $tmp/db 
    valgrind $synder \
        -i "$PWD/test/test-data/multi-block/c.gff" \
        -s $tmp/db/a_b.txt \
        -c search > /dev/null 2> valgrind.log
    valgrind_exit_status=$?
    rm -rf tmp
fi


#---------------------------------------------------------------------

total=$(( total_passed + total_failed))
emphasize "$total_passed tests successful out of $total"

if [[ $valgrind_checked == 0 ]]
then
    warn "valgrind not found, no memory tests performed\n"
else
    if [[ $valgrind_exit_status == 0 ]]
    then
        emphasize_n "valgrind pure"
        echo " (for synder search of multi-block/c.gff against multi-block/map.syn)"
        rm valgrind.log
    else
        warn "valgrind failed - see valgrind.log\n"
    fi
fi

if [[ $total_failed > 0 ]]
then
    warn "$total_failed tests failed\n"
    exit 1
else
    exit 0
fi
