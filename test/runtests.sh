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
    echo -e $o
}

emphasize(){
   emphasize_n "$1" 
   echo
}

emphasize_n(){
    [[ -t 1 ]] && o="\e[1;39m$1\e[0m" || o=$1
    echo -ne $o
}

runtest(){
    dif=$1
    base=$2
    msg=$3
    tmp=/tmp/synder-$RANDOM$RANDOM
    mkdir $tmp
    echo -n "Testing $msg ... "
    $synder -d $dir/map.syn a b $tmp/db 
    $synder -i $dir/$base.gff -s $tmp/db/a_b.txt -c search > $tmp/a
    diff $tmp/a $dir/${base}-exp.txt > /dev/null 
    if [[ $? == 0 ]]
    then
        total_passed=$(( $total_passed + 1 ))
        echo "OK"
    else
        warn "FAIL"
        total_failed=$(( $total_failed + 1 ))
        echo "======================================="
        emphasize "expected output:"
        cat $dir/${base}-exp.txt
        emphasize "observed output:"
        cat $tmp/a
        emphasize "query gff:"
        cat $dir/$base.gff
        emphasize "synteny map:"
        cat $dir/map.syn
        echo -e "---------------------------------------\n"
    fi

    rm -rf $tmp
}

#---------------------------------------------------------------------
dir="$PWD/test/test-data/one-block"
announce "\nTesting with synteny map length == 1"
runtest $dir hi     "query downstream of block"
runtest $dir within "query within block"
runtest $dir lo     "query upstream of block"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/two-block"
announce "\nTesting with synteny map length == 2"
runtest $dir hi      "query downstream of all blocks"
runtest $dir between "query between two blocks"
runtest $dir lo      "query upstream of all blocks"

#---------------------------------------------------------------------
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
    warn "valgrind not found, no memory tests performed"
else
    if [[ $valgrind_exit_status == 0 ]]
    then
        emphasize_n "valgrind pure"
        echo " (for synder search of multi-block/c.gff against multi-block/map.syn)"
        rm valgrind.log
    else
        warn "valgrind failed - see valgrind.log"
    fi
fi

if [[ $total_failed > 0 ]]
then
    warn "$total_failed tests failed"
    exit 1
else
    exit 0
fi
