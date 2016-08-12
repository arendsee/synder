#!/bin/bash
set -u
set -e

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
  -z input is 1-based (default, 0-based)
OPTIONAL ARGUMENTS
  -q query chromosome length file
  -t target chromosome length file
INPUT FORMAT
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

base=0
query= target= outdir= input= query_gf= target_gf=
while getopts "ha:b:d:i:t:q:z" opt; do
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
        q)
            query_gf="$OPTARG" 
            if [[ ! -f "$query_gf" ]]
            then
                die "Cannot read query genome length file '$query_gf'"
            fi
            ;;
        t)
            target_gf="$OPTARG" 
            if [[ ! -f "$target_gf" ]]
            then
                die "Cannot read target genome length file '$target_gf'"
            fi
            ;;
        d)
            outdir="$OPTARG"
            [[ -d $outdir ]] || mkdir $outdir || die "Cannot access directory '$outdir'"
            ;;
        z)
            base=1
            ;;
        ?)
            die "Unrecognized argument"
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
# Parse utilities
# ============================================================================

# Add two columns to the input table
#  1. Map of sequence names to 0-based indices (alphabetically)
#  2. Map of intervals on a sequence to 0-based indices (by start and stop) 
# Arguments:
# ARG1: sequence id (e.g. Chr1 or scaffold_1) 
# ARG2: start position relative to sequence id $1
# ARG3: stop position relative to sequence id $1
append_counts() {
sort -k${1},${1} -k${2},${2}n -k${3},${3}n | awk -v NC=$1 '
        BEGIN { blkid=0; seqid=0 }
        NR == 1 { seqname = $NC }
        seqname != $NC {
            seqid++
            blkid = 0
            seqname = $NC
        }
        { 
            print $0, seqid, blkid 
            blkid++
        }
    '
}


# Map sequence ids to names
#
# Lines beginning in '>' indicate the species name and contig count. In the
# example below, the name species is Arabidopsis_thaliana and it has 7
# chromosomes.
#
# The lines beginning with '$' contain four columns:
#  1) contig id
#  2) number of intervals on the contig
#  3) name of the contig
#  4) contig length
#
# An '@' sign terminates the species entry
#
# Example:
# >  Arabidopsis_thaliana  7
# $  0  60721  Chr1   30427671
# $  1  37435  Chr2   19698289
# $  2  44834  Chr3   23459830
# $  3  33059  Chr4   18585056
# $  4  53222  Chr5   26975502
# $  5      8  ChrM     366924
# $  6    283  ChrC     154478
# @

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

add_chromosome_length() {
    awk '
        BEGIN{ FS="\t"; OFS="\t" }
        NR == FNR { a[$1] = $2; next }
        $1 == "$" {
            if( $4 in a ) {
                $0 = $0 "\t" a[$4]
            } else {
                errmsg = $4 " from synteny file not in scaffold length file" 
                exit 1
            }
        }
        { print }
        END{
            if(errmsg){
                print errmsg > "/dev/stderr" 
            }
        }
    ' $1 /dev/stdin
}

add_default_chromosome_length() {
    awk '
        $1 == "$" { $0 = $0 "\t" 1000000000 }
        { print }
    ' /dev/stdin
}

parse() {
    input=$1   # synteny map
    outbase=$2 # output file basename
    query=$3   # string: query name (e.g. Arabidopsis_thaliana)
    target=$4  # string: target name
    outdb=${outbase}.txt       # ultimate output
    outtmp=${outbase}_temp.txt # temporary output
    # outtmp columns (1-8 are from the synteny map, see USAGE):
    # 10 qseqid - index for query sequence
    # 11 qblkid - index for interval in a query sequence array
    # 12 tseqid - index for target sequence
    # 13 tblkid - index for interval in a target sequence array
    cat $input |
        # Append qseqid and qblkid
        append_counts 1 2 3 | 
        # Append tseqid and tblkid
        append_counts 4 5 6 > $outtmp

    if [[ -r $query_gf ]]
    then
        write_side $query 1 9 < $outtmp |
            add_chromosome_length $query_gf > $outdb
    else
        write_side $query 1 9 < $outtmp |
            add_default_chromosome_length > $outdb
    fi

    if [[ -r $target_gf ]]
    then
        write_side $target 4 11 < $outtmp |
            add_chromosome_length $target_gf >> $outdb
    else
        write_side $target 4 11 < $outtmp |
            add_default_chromosome_length >> $outdb
    fi


    # Output link data and perform offsets
    awk -v base=$base '
        {
            qseqid = $9
            qblkid = $10
            qstart = $2   - base
            qstop  = $3   - base
            tseqid = $11
            tblkid = $12
            tstart = $5   - base
            tstop  = $6   - base
            strand = $8
            print "$", qseqid, qblkid, qstart, qstop, tseqid, tblkid, tstart, tstop, strand
        }
    ' $outtmp |
        sort -k3,3n -k4,4n >> $outdb
    rm $outtmp
}


# ============================================================================
# Now run the thing
# ============================================================================

parse $input "${outdir}/${query}_${target}" $query $target
