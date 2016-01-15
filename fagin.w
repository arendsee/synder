\documentclass{cweb}
\usepackage{rcs}
\usepackage{enumerate}

\def\fagin{{\tt Fagin\/}}

\begin{document}


\title{Fagin}
\author{Zebulun Arendsee}
\RCSdate $Date: 1995/08/29 17:27:57 $

\maketitle

% not very interesting, only one starred section
% \tableofcontents

@* Introduction to \fagin{}.

Fagin is a tool designed to integrate sequence match, syntenic, and
transcriptomic data to trace the history of genes having obscure or unknown
lineage (ghouls).

Traditionally, BLAST is used to identify orphan genes. However, orphans so
identified are a heterogenous group, containing ancient, rapidly evolving
genes; true de novo genes; small genes genes BLAST simply misses; genes that
are not annotated as genes in other species. Some papers perform this task with
a pipeline of tools mixed with manual review. My goal is to build a formal,
broadly applicable program.

@*Preamble.

@c

#include <stdio.h>
#include <stdlib.h>

@*Input.

@*2 GFF Input.

A General Feature Format (GFF) file contains the locations of features relative
to a string, usually a biological sequence.

@c

#define SEQID_BUFFER_LENGTH 32
#define ID_BUFFER_LENGTH 32
#define TYPE_BUFFER_LENGTH 32

typedef struct GFFEntry
{
    char seqid[SEQID_BUFFER_LENGTH];
    char type[TYPE_BUFFER_LENGTH];
    int start;
    int stop;
    char strand;
    char id[ID_BUFFER_LENGTH];
} GFFEntry;

void load_gff(FILE *fp, GFFEntry * gff){
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    fseek(fp, 0, SEEK_SET);

    int i = 0;

    while ((read = getline(&line, &len, fp)) != -1) {
        if(line[0] == '#')
            continue;
        sscanf(line,
               "%s %*s %s %d %d %*c %c %*c %s",
               &gff[i].seqid, &gff[i].type, &gff[i].start, &gff[i].stop, &gff[i].strand, &gff[i].id);
        i++;
    }
}

@*2Synteny Input.

Synteny files match intervals in one string to intervals in another string.
They also contain a score for the match (percent identity) and specifiy the
direction/sense of the match (strand).

They should have the following columns in exactly the following order:

\begin{enumerate}
    \item  query id [string]

    \item  query start [int]

    \item  query stop [int]

    \item  tardet id [string]

    \item  target stop [int]

    \item  target start [int]

    \item  percent identity [float]

    \item  strand [char], this can be '+', '-', or '.' (if unknown/irrelevant)
\end{enumerate}

@c

typedef struct SynEntry
{
    char *qseqid;
    char *tseqid;
    int qstart;
    int qend;
    int tstart;
    int tend;
    float pident;
    char strand;
} SynEntry;

@*Utilities.

These are a set of functions that do obvious things.

@*2count\_lines.

Count the lines in a file given a file handle.

@c
int count_lines(FILE *fp){
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int nlines = 0;

    int initial = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    while ((read = getline(&line, &len, fp)) != -1) {
        if(line[0] != '#')
            nlines++;
    }

    fseek(fp, initial, SEEK_SET);

    return nlines;
}

@*Main Function.

And the big boy.

@c

int main(int argc, char* argv[])
{
    FILE *synf;    // synteny file
    int nsyn = 0;  // lines in synteny file

    FILE *qgfff;   // query gff file
    int ngff = 0;  // lines in gff file

    int i = 0; // generic index

    // There must be 2 arguments:
    // ARG1: synteny filename
    // ARG2: gff filename
    if(argc != 3){
        return 1;
    }

    synf = fopen(argv[1], "rb");
    qgfff = fopen(argv[2], "rb");

    nsyn = count_lines(synf);
    /* struct SynEntry syn[nsyn]; */

    ngff = count_lines(qgfff);
    GFFEntry * gff = (GFFEntry *)malloc(ngff * sizeof(GFFEntry));

    load_gff(qgfff, gff); 

    for(i = 0; i < 10; i++){
        printf("%s\t%s\t%d\t%d\t%c\t%s\n",
               &gff[i].seqid, &gff[i].type, gff[i].start, gff[i].stop, gff[i].strand, &gff[i].id);
    }

    fclose(synf);
    fclose(qgfff);

    free(gff);

    return 0;
}

@

\bibliographystyle{plain}
\bibliography{fagin}

@

\cwebIndexIntro{%
    Here is a list of the identifiers used, and where they appear.
Underlined entries indicate the place of definition. Error messages
are also shown.
    }

\end{document}
