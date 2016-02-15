#!/bin/bash
set -u


# ============================================================================
# Parse Arguments
# ============================================================================

print_help(){
    cat << EOF
prepare synteny file for fagin
REQUIRED ARGUMENTS
  -a name of query organism
  -b name of target organism
  -i name of input sequence
  -d output data directory
INPUT FORMAT (read from STDIN)
  qseqid qstart qend tseqid tstart tend score strand
Where:
  qseqid: name of the query sequence
  qstart: start position of the query
  qstop:  stop position of the query
  tseqid: name of the target sequence
  tstart: start position of the target
  tstop:  stop position of the target
  score:  score of the alignment (high is good)
  strand: '+' for same strand, '-' for opposing strand 
EOF
}

die(){
    echo "ERROR: $1" >&2
    exit 1
}

die-instructively(){
    echo -e '\e[1;31mERROR: '$1'\e[0m' >&2
    print_help >&2
    exit 1
}



# ============================================================================
# Parse utilities
# ============================================================================

append_counts() {
sort -k${1},${1}n -k${2}n | awk -v NC=$1 '
        BEGIN { blkid=0; seqid=0 }
        NR == 1 { s = $NC }
        s != $NC {
            seqid++
            blkid = 0
            s = $NC
        }
        { 
            print $0, seqid, blkid 
            blkid++
        }
    '
}

write_side() {
    awk -v seqname=$1 -v namecol=$2 -v idcol=$3 '
        BEGIN{ OFS="\t" }
        {
            counts[$idcol]++
            map[$idcol] = $namecol
        }
        END{
            for (k in map){ nseqs++ }
            print ">", seqname,  nseqs
            for (k in map){
                print "$", k, counts[k], map[k]
            }
            print "@"
        }
    '
}

parse() {
    intput=$1
    outbase=$2
    query=$3
    target=$4
    outdb=${outbase}.txt
    outtmp=${outbase}_temp.txt
    append_counts 1 2 < $input |
        append_counts 4 5 |
        awk '{print $0, linkid++}' > $outtmp
    > $outdb
    write_side $query  1 9  < $outtmp >> $outdb
    write_side $target 4 11 < $outtmp >> $outdb
    awk '{print "$", $9, $10, $2, $3, $11, $12, $5, $6, $13}' $outtmp |
        sort -k1,1n -k2,2n >> $outdb
    rm $outtmp
}


# ============================================================================
# Parse commandline 
# ============================================================================

query= target= outdir= input=
while getopts "ha:b:d:i:" opt; do
    case $opt in
        h)
            print_help
            exit 0 ;;
        a)
            query="$OPTARG" ;;
        b)
            target="$OPTARG" ;;
        i)
            input="$OPTARG" 
            if [[ ! -f $input ]]
            then
                die "Cannot read input file '$input'"
            fi
            ;;
        d)
            outdir="$OPTARG"
            [[ -d $outdir ]] || mkdir $outdir || die "Cannot access directory '$outdir'"
            ;;
    esac
done
# Check for existance of required arguments
for arg in query target outdir input; do
    if [[ -z ${!arg} ]]; then
        die-instructively 'Missing required argument'
    fi
done



# ============================================================================
# Now run the thing
# ============================================================================

parse $input "${outdir}/${query}_${target}" $query $target
