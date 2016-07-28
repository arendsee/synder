#!/usr/bin/env bash
set -u

synder=$PWD/synder

total_passed=0
total_failed=0

warn(){
    echo -e "\e[1;31m$1\e[0m"
}

emphasize(){
   emphasize_n "$1" 
   echo
}

emphasize_n(){
    echo -ne "\e[01;39m$1\e[0m"
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
dir="$PWD/test/A"
runtest $dir hi     "itree (set A - hi)"
runtest $dir within "itree (set A - within)"
runtest $dir lo     "itree (set A - lo)"

#---------------------------------------------------------------------
dir="$PWD/test/B"
runtest $dir hi "itree (set B - hi)"
runtest $dir within "itree (set B - within)"
runtest $dir lo "itree (set B - lo)"

#---------------------------------------------------------------------
dir="$PWD/test/C"
runtest $dir a "itree (set C - a)"
runtest $dir b "itree (set C - b)"
runtest $dir c "itree (set C - c)"
runtest $dir d "itree (set C - d)"
runtest $dir e "itree (set C - e)"
runtest $dir f "itree (set C - f)"
runtest $dir g "itree (set C - g)"
runtest $dir h "itree (set C - h)"
runtest $dir i "itree (set C - i)"
runtest $dir j "itree (set C - j)"

#---------------------------------------------------------------------
echo
 
#---------------------------------------------------------------------
# valgrind

valgrind_checked=0
valgrind_exit_status=1
if hash valgrind 2> /dev/null; then
    valgrind_checked=1
    dir="$PWD/test/C"
    tmp=/tmp/synder-$RANDOM
    mkdir $tmp
    $synder -d $dir/map.syn a b $tmp/db 
    valgrind $synder \
        -i "$PWD/test/C/c.gff" \
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
        echo " (for synder search of C/c.gff against C/map.syn)"
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
