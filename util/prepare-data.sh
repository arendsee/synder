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
# Parse data into four files 
#  * <query>_<target>_left.tab   - the query, or input
#  * <query>_<target>_right.tab  - the target
#  * <query>_<target>_score.tab  - a list of scores linked to blockid
#  * <query>_<target>_strand.tab - a list of the direction of the matches
# ============================================================================
parse(){
    input=$1
    outbase=$2
    sort $input -k1 -nk2,2 | awk -v OFS="\t" -v base=$outbase '
        {
            qseqid = $1
            qstart = $2
            qstop  = $3
            tseqid = $4
            tstart = $5
            tstop  = $6
            score  = $7
            strand = $8

            blockid++
            qcounts[qseqid]++
            lines[blockid] = qstart "\t" qstop

            if(NR == 1 || $1 != prev){
                start[blockid] = qseqid
            }
            prev = $1

            # Record block data
            print blockid, score > base "_score.tab"
            print blockid, strand > base "_strand.tab"

            # send target data to STDOUT
            print blockid, tseqid, tstart, tstop
        }
        END {
            out = base "_left.tab"
            for(i = 1; i <= blockid; i++){
                if(i in start){
                    print "> " ++j " " start[i] " " qcounts[start[i]] > out
                }
                print i, lines[i] > out
            }
        }
    ' | sort -k1 -nk3,3 | awk -v OFS="\t" -v base=$outbase '
        {
            blockid = $1
            seqid = $2
            start = $3
            stop = $4

            lines[blockid] = start "\t" stop
            qcounts[seqid]++
            if(NR == 1 || seqid != prev){
                seq_start[blockid] = seqid
            }
            prev = seqid
        }
        END {
            out = base "_right.tab"
            for(k in lines){
                if(k in seq_start){
                    print "> " ++j " " seq_start[k] " " qcounts[seq_start[k]] > out
                }
                print k, lines[k] > out
            }
        }
    '
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
            if [[ ! -d $outdir ]]
            then
                die "Cannot access directory '$outdir'"
            fi
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

parse $input "${outdir}/${query}_${target}"
