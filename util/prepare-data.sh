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
# The lines beginning with '$' contain three columns: 1) the contig id, 2) the
# number of intervals on the contig, and 3) the name of the contig.
#
# An '@' sign terminates the species entry
#
# Example:
# >	Arabidopsis_thaliana	7
# $	0	283	chloroplast
# $	1	60721	Chr1
# $	2	37435	Chr2
# $	3	44834	Chr3
# $	4	33059	Chr4
# $	5	53222	Chr5
# $	6	8	mitochondrion
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
    # 14 linkid - index for the link between a query and target interval
    cat $input |
        # Append qseqid and qblkid
        append_counts 1 2 3 | 
        # Append tseqid and tblkid
        append_counts 4 5 6 | 
        # Append linkid
        awk '{print $0, linkid++}' > $outtmp
    write_side $query 1 9 < $outtmp >  $outdb
    write_side $target 4 11 < $outtmp >> $outdb
    awk '
        {
            qseqid=$9
            qblkid=$10
            qstart=$2
            qstop=$3
            tseqid=$11
            tblkid=$12
            tstart=$5
            tstop=$6
            linkid=$13
            strand=$8
            print "$", qseqid, qblkid, qstart, qstop, tseqid, tblkid, tstart, tstop, linkid, strand
        }
    ' $outtmp |
        sort -k3,3n -k4,4n >> $outdb
    rm $outtmp
}


# ============================================================================
# Now run the thing
# ============================================================================

parse $input "${outdir}/${query}_${target}" $query $target
