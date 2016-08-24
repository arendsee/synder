#!/usr/bin/env bash
set -u

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

usage (){
cat << EOF
Test the final output of synder
OPTIONAL ARGUMENTS
  -h  print this help message
  -x  die on first failure 
  -d  print full debugging info and link input file (g), database file (d),
      observed (o) and expected (e) files to the working directory.
  -v  verbose
  -o  redirect gdb output to this file (or tty)
  -m  test memory with valgrind
EOF
    exit 0
}

die_on_failure=0
debug=0
valgrind=0
verbose=0
gdb_out="none"
while getopts "hdxvmo:" opt; do
    case $opt in
        h)
            usage ;;
        d) 
            debug=1 ;;
        x)
            die_on_failure=1 ;;
        v)
            verbose=1 ;;
        o)
            gdb_out=$OPTARG
            ;;
        m)
            if [[ ! -z `type -P valgrind` ]]
            then
                valgrind=1
            else
                warn "Valgrind not found"   
            fi
            ;;
        ?)
            echo "Illegal arguments"
            exit 1
    esac 
done

synder=$PWD/synder

total_passed=0
total_failed=0
valgrind_checked=0
valgrind_exit_status=0

# A function to select which parts of the output should be compared
# Since flags are currently in flux, test only the first 7 columns
filter () {
    sort | cut -f1-8,10-12
}

filter_plus_one () {
    awk -v OFS="\t" '{$3++ ; $4++ ; $6++ ; $7++ ; print}' | filter
}

# This variable will be set for each set of test
# It specifies where the input files can be found
dir=
runtest(){
    base=$1
    msg=$2
    errmsg=${3:-0}
    out_base=${4:-0}

    # Write output to this folder, if all goes well, it will be deleted
    tmp=/tmp/synder-$RANDOM$RANDOM
    mkdir -p $tmp/db

    # initialize temporary files
    gff=$dir/$base.gff
    map=$dir/map.syn
    val=$tmp/v
    cmd=$tmp/c
    obs=$tmp/o
    exp=$tmp/e
    tdb=$tmp/db/a_b.txt

    echo -n "Testing $msg ... "

    # Default synder database building command
    db_cmd="$synder -d $map a b $tmp/db"

    # Query genome length file
    tgen=$dir/tgen.tab
    # Target genome length file
    qgen=$dir/qgen.tab

    # If the length files are given, add to db command
    if [[ -f $tgen ]]; then
        db_cmd="$db_cmd $tgen"
    fi
    if [[ -f $tgen && -f $qgen ]]; then
        db_cmd="$db_cmd $qgen"
    fi

    if [[ $out_base == 1 ]]
    then
        cat $dir/${base}-exp.txt | filter_plus_one > $exp
    else
        cat $dir/${base}-exp.txt | filter > $exp
    fi

    # Build database
    $db_cmd

    if [[ $? != 0 ]]
    then
        total_failed=$(( $total_failed + 1 ))
        warn "FAILED - Could not build database\n"
        echo "failed command:"
        echo $db_cmd
    else

        if [[ $debug -eq 1 ]]
        then
            ln -sf $gff  g 
            ln -sf $map  m
            ln -sf $exp  e
            ln -sf $obs  o
            ln -sf $tdb  d
            ln -sf $val  v
            ln -sf $cmd  c
        fi

        synder_cmd=$synder

        [[ $out_base == 1 ]] && synder_cmd="$synder_cmd -b 0011 "
        synder_cmd="$synder_cmd -i g"
        synder_cmd="$synder_cmd -s d"
        synder_cmd="$synder_cmd -c search"

        # command for loading into gdb
        echo "set args $synder_cmd"  >  $cmd
        echo "source .gdb_cmds"     >> $cmd
        if [[ $gdb_out != "none" ]]
        then
            echo "set logging off"           >> $cmd
            echo "set logging file $gdb_out" >> $cmd
            echo "set logging redirect on"   >> $cmd
            echo "set logging on"            >> $cmd
            echo "gdb -tui --command c" > x
            chmod 755 x
        fi

        # Ensure all input files are readable
        for f in $gff $exp $map $tdb;
        do
            if [[ ! -r "$f" ]]
            then
                warn "Input file $f not readable"
            fi
        done

        $synder_cmd 2>&1 > $obs
        if [[ $? != 0 ]]
        then
            warn "Synder terminated with non-zero status:\n"
            echo $synder_cmd
            exit 1
        fi

        $synder_cmd | filter > $obs

        if [[ $valgrind -eq 1 ]]
        then
            valgrind $synder_cmd > /dev/null 2> $val
            valgrind_exit_status=$?
        fi

        diff $obs $exp > /dev/null 

        if [[ $? == 0 ]]
        then
            if [[ $valgrind -eq 1  &&  $valgrind_exit_status -ne 0 ]]
            then
                warn "memory bug"
                [[ $die_on_failure -eq 1 ]] && exit 1
            else
                # clear all temporary files
                rm -rf x g e o d v /tmp/synder-*
                total_passed=$(( $total_passed + 1 ))
                echo "OK"
            fi
        else
            warn "FAIL"
            echo " (in `basename $dir`/)"
            [[ $errmsg == 0 ]] || (echo -e $errmsg | fmt)
            total_failed=$(( $total_failed + 1 ))
            echo "======================================="
            emphasize_n "expected output"; echo ": (${base}-exp.txt)"
            cat $exp
            emphasize "observed output:"
            cat $obs
            emphasize_n "query gff"; echo ": (${base}.gff)"
            column -t $gff
            emphasize_n "synteny map"; echo ": (map.syn)"
            column -t $map
            if [[ $debug -eq 1 && $verbose -eq 1 ]]
            then
                echo "Debugging files:"
                echo " * g - input GFF file"
                echo " * o - observed output"
                echo " * e - expected output"
                echo " * d - database"
                echo " * v - valgrind report (if valgrind is present)"
                echo "Synder database command:"
                echo $db_cmd
                echo "Synder command:"
                echo $synder_cmd
            fi
            echo -e "---------------------------------------\n"

            [[ $die_on_failure -eq 1 ]] && exit 1

        fi
    fi
}

#---------------------------------------------------------------------
dir="$PWD/test/test-data/one-block"
announce "\nTesting with synteny map length == 1"
runtest hi     "query after of block"
runtest within "query within block"
runtest lo     "query before of block"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/two-block"
announce "\nTesting with synteny map length == 2"
runtest hi      "query downstream of all blocks"
runtest between "query between two blocks"
runtest lo      "query upstream of all blocks"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/multi-block"
announce "\nTesting with 5 adjacent blocks on the same strand"
runtest a "extreme left"
runtest b "inbetween two adjacent blocks"
runtest c "starts inbetween adjacent blocks"
runtest d "stops inbetween adjacent blocks"
runtest e "inbetween two adjacent blocks"
runtest f "starts before block 3, ends after block 3"
runtest g "starts in block 2, ends after block 3"
runtest h "starts before block 2, ends after block 3"
runtest i "starts in block 2, ends in block 2"
runtest j "extreme right"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/simple-duplication"
announce "\nTest simple tandem duplication"
runtest between "query starts between the duplicated intervals"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/one-interval-inversion"
announce "\nTest when a single interval is inverted"
runtest between "query next to inverted interval"
runtest over    "query overlaps inverted interval"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/two-interval-inversion"
announce "\nTest when two interval are inverted"
runtest beside   "query next to inverted interval"
runtest within   "query between inverted intervals"
runtest spanning "query spans inverted intervals"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/tiny-indel-query-side"
announce "\nTest when a small interval interupts on one side"
runtest beside "query side"
dir="$PWD/test/test-data/tiny-indel-target-side"
runtest beside "target side"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/tandem-transposition"
announce "\nTest tandem transposition"
runtest beside "query beside the transposed pair"
runtest within "query between the transposed pair"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/irregular-overlaps"
announce "\nTest target side internal overlaps"
runtest left "left side" "You are either 1) not sorting the by_stop vector
in Contig by Block stop positions, or 2) are snapping the search interval left
boundary to a Block that is nearest by start, but not be stop."
runtest right "right side"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/multi-chromosome"
announce "\nTest two intervals on same query chr but different target chr"
runtest between "between the query intervals"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/inverted-extremes"
announce "\nExtreme value resulting from an inversion"
runtest extreme "between the query intervals, extreme SI"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/deletion"
announce "\nDeletion tests (adjacent bounds in target)"
runtest between "query is inbetween"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/off-by-one"
announce "\nTest overlap edge cases"
runtest a "overlap of 1"

# #---------------------------------------------------------------------
# # TODO Find a good way to deal with this case:
# dir="$PWD/test/test-data/synmap-overlaps"
# announce "\nsyntenic overlaps"
# runtest simple "Between the weird"

#---------------------------------------------------------------------
dir="$PWD/test/test-data/unassembled"
announce "\nMappings beyond the edges of target scaffold"
runtest lo "query is below scaffold"
runtest adj-lo "query is just below the scaffold"
runtest adj-hi "query is just above the scaffold"
runtest hi "query is above the scaffold"

runtest lo "test with 1-base" 0 1

#---------------------------------------------------------------------
echo

total=$(( total_passed + total_failed))
emphasize "$total_passed tests successful out of $total"

if [[ $total_failed > 0 ]]
then
    warn "$total_failed tests failed\n"
    exit 1
else
    exit 0
fi
