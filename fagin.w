\documentclass{cweb}
\usepackage{enumerate}
\def\fagin{{\tt Fagin\/}}

\begin{document}
\title{Fagin}
\author{Zebulun Arendsee}
\maketitle
\tableofcontents

@*Introduction to \fagin{}. 

\fagin{} is a tool designed to integrate sequence match, syntenic, and
transcriptomic data to trace the history of genes having obscure or unknown
lineage (ghouls).

Traditionally, BLAST is used to identify orphan genes. However, orphans so
identified are a heterogenous group, containing ancient, rapidly evolving
genes; true de novo genes; small genes genes BLAST simply misses; genes that
are not annotated as genes in other species. Some papers perform this task with
a pipeline of tools mixed with manual review. My goal is to build a formal,
broadly applicable program.

Traditionally, in papers such as xxx de novo orphans have been identified by
passing all genes through a series of filters. Usually something like this:

\begin{enumerate}

    \item All proteins are blasted against some protein database, those with no
    extraspecies homologs are classified as "orphans".

    \item The DNA coding sequence of each orphan is blasted against the full
    genomes of several relatives of the focal species. Orphans with hits to
    things that could be valid ORFs are discarded.

    \item Other filters (skim lit for this) ...

\end{enumerate}

However there are several limitations of this approach.

\begin{itemize}

    \item The process is not encapsulated in a broadly applicable pipeline, so
    it taxes researcher time.

    \item It is reliant of DNA sequence similarity of (usually) very short
    genes. Homologs may not be identified even if they do exist.

    \item It identifies a small, conservative set of de novo orphans and
    discards the rest as unknown.

\end{itemize}

\fagin{} is designed to sort all ghouls into distinct classes. Homologous
genomic intervals can be identified in relatives by rooting the ghoul region to
its context using genome wide synteny predictions (from Satsuma for example).

<insert inkscape diagram showing mappings>

Genes are classified into origin stories based on the following features:
\begin{description}
    \item[overlaps syntenic block] The gene is anchored to an exact position.
    \item[context is defined] 
    \item[context is missing] The gene matches a section of the genome that may not have been sequenced.
    \item[syntenic match to a gene]
    \item[syntenic match to a stranscript CDS]
    \item[syntenic match to a genomic ORF]
\end{description}

<insert more inkscape diagrams>

Here are a few possible origin stories:
\begin{description}
    \item[de novo orphan] The gene maps to a non-coding region
    \item[fast evolver] The gene maps to a possible gene, but one that evolves so quickly similarity is lost.
    \item[uncertain] The synteny is convoluted, classification failed
    \item[missing data] Context includes unknown, can't infer
    \item[false orphan] Annotationn missing in cousin
\end{description}

I will not always be able to infer class with certainty. Will need to develop a
statistical model to estiate uncertainty in classification.

@*1Program overview. 

Here is a stub main function

@c
@<Include statements@>@/
@<Global declarations@>@/
@<Function prototypes@>@/
@<The main program@>
@<Data structure functions@>@/
@<Program input functions@>@/
@



% =============================================================================
@*1Main function.

Currently this system merely loads a tree and prints it.

@<The main program@>=

int main(int argc, char* argv[]){

    if(argc < 2){
        printf("USAGE: ./a.out <filename>\n");
        exit(EXIT_FAILURE);
    }

    // get doubly-linked list of nodes
    Node * root = loadNodeList(argv[1]);

    // raise this list into a tree
    buildTree(root, NLEVELS);

    // wire siblings
    wireGenomeTree(root);

    recursivePrintNode(root);

    freeNode(root);

    exit(EXIT_SUCCESS);
}
@



% =============================================================================
@*1Dependencies.

We must include the standard I/O definitions, since we want to send formatted
output to |stdout| and |stderr|.

@<Include statements@>=
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

@

@<Global declarations@>=
#define NAME_LENGTH 20

#define GENOME 0
#define CHR    1
#define GENE   2
#define MRNA   3
#define CDS    4

#define NLEVELS 5

typedef struct Node {
    char name[NAME_LENGTH];
    unsigned char type;
    size_t size;
    struct Node ** children;
    struct Node * next;
    struct Node * last;
} Node;
@



% =============================================================================
@*1Data structure.

@<Data structure functions@>=
  @<Node print functions@>@/
  @<Node memory functions@>@/
  @<Node build from file@>@/
@

% -----------------------------------------------------------------------------
@*2 Node declaration and default values.


% -----------------------------------------------------------------------------
@*2 Node print functions.

@<Node print functions@>=
void printNode(Node * node){
    printf("%s\t\ttype=%d size=%d last=%s next=%s\n",
            node->name,
            node->type,
            node->size,
            node->last->name,
            node->next->name);
}

void recursivePrintNode(Node * node){
    printNode(node);
    for(unsigned int i = 0; i < node->size; i++){
        if(node->children != NULL)
            recursivePrintNode(node->children[i]);
    }
}
@


% -----------------------------------------------------------------------------
@*2 Node memory allocation.

@<Node memory functions@>=

// Allocate space for a new node and set defaults
Node * newNode(){
    Node * nd = (Node *) malloc(sizeof(Node));
    nd->type = 99;
    nd->size = 0;
    nd->children = NULL;
    nd->next = NULL;
    nd->last = NULL;
    return(nd);
}

void freeNode(Node * node){
    for(int i = 0; i < node->size; i++){
        freeNode(node->children[i]);
    }
    free(node->children);
    free(node);
}

void freeList(Node * node){
    if(node->next != NULL){
        freeList(node->next);
    }
    freeNode(node);
}
@

% -----------------------------------------------------------------------------
@*1 Build a Genomic Christmas tree.

This function builds a tree from a doubly-linked list of Node structures. The
nodes are assumed to be in a depth-first recursive order.
     
The output data structure is a tree augmented with doubly-linked lists joining
certain branches. It can be searched recursively from root (genome) to leaf
(exon). Two levels, chromosome and gene, can also be iterated through with LAST
and NEXT links. This allows direct access to siblings.  It also permits arrays
of pointers to nodes to be easily copied for analysis, sorting, and whatnot
outside the main data structure. Genes do not link across chromosomes, which
prevents false inferences of adjacency.

Building the tree requires three passes.

@<Node build from file@>=
Node * buildTree(Node * root, size_t nlevels){
    @<buildTree: local variables@>@/
    @<buildTree: count children@>@/
    @<buildTree: assign children to parent arrays@>@/
    return(root);
}
@<Node internal doubly-linked listing@>@/
@

@<buildTree: local variables@>=
    Node * node = root->next;

    size_t counts[nlevels];
    Node * parents[nlevels];
    for(int i = 0; i < nlevels; i++){
        counts[i] = 0;
        parents[i] = NULL;
    }
    parents[0] = root;
@

@*2 First pass: count the children.

The first pass over the list simply counts the number of children descending
from each node.

@<buildTree: count children@>=
    while(node != NULL){
        parents[node->type - 1]->size++;
        parents[node->type] = node;
        node = node->next; 
    }
@

@*2 Second pass: allocate arrays and assign children.

The second pass allocates a pointer array of the size calculated in the first
pass to each node. It fills each array with pointers to each of the nodes
children.

@<buildTree: assign children to parent arrays@>=
    node = root;
    size_t pid;
    while(node != NULL){

        // Parent's level (e.g. 0 is GENOME and 4 is EXON)
        pid = node->type - 1;

        // make container for children
        if(node->size > 0)
            node->children = (Node**)malloc(node->size * sizeof(Node*));

        // start collecting children for this node at the appropriate level
        parents[node->type] = node;
        // reset child count at this level
        counts[node->type] = 0;

        // skip root
        if(node->last != NULL){
            // add this child to its parent's children array
            parents[pid]->children[counts[pid]] = node;

            // count this child (needed for correct assignment of next sibling)
            counts[pid]++;  
        }

        // progress to next node (remember this is still a flat, doubly-linked list)
        node = node->next;
    }
@

@*2 Third pass: building stratified linked-lists within the tree.

This final pass is optional. The method chosen will vary depending on the
category of tree.

Determine the last/next relationships for a tree of genome data

This function assumes the GENOME $\rightarrow$ CHR $\rightarrow$ GENE
$\rightarrow$ MRNA $\rightarrow$ CDS/exon structure. The actual type names are
not important, since they are mapped to numbers (GENOME == 0, CDS == 4).

@<Node internal doubly-linked listing@>=
void wireGenomeTree(Node * node){
    for(unsigned int i = 0; i < node->size; i++){
        wireGenomeTree(node->children[i]);
    }
    if(node->size > 0){
        node->children[0]->last = NULL;
        node->children[node->size-1]->next = NULL;
    }
    if(node->size > 1 && (node->type == GENOME || node->type == CHR)){
        for(unsigned int i = 1; i < node->size; i++){
            node->children[i]->last = node->children[i-1];
            node->children[i-1]->next = node->children[i];
        }
    }
}
@


@*1 Program input functions.

@*2 Loading genome features.

A General Feature Format (GFF) file contains the locations of features relative
to a string, usually a biological sequence.

\fagin{} does not directly read GFF files, however. It requires a very specific
input that exactly follows a specific order. It must have at least three space
separated arguments: 

For example:

\begin{verbatim}
chr     1      Chr1     1   10000
gene    2        g1  3631    5899
mRNA    3      g1m1  3631    5899
CDS     4    g1m1c1  3760    3913
CDS     4    g1m1c2  3996    4276
CDS     4    g1m1c3  4486    4605
gene    2        g2  5928    8737
...
\end{verbatim}

\begin{enumerate}

  \item TYPE This is actually ignored, but must be present.

  \item LEVEL \fagin{} uses this number to determine hierarchy level. 0 is root,
  but mmust not be included.

  \item NAME A label with no spaces, it must be less than NAME\_LENGTH
  characters long. There are no other constraints (currently).  It needn't be
  unique and can include any non-whitespace ASCII character.

  \item BEG The starting position of the interval.

  \item END The ending position of the interval.
\end{enumerate}

The file may contain additional columns of data, but at present these will
be ignored.

The order is REQUIRED to be depth-first recursion order based on the TYPE
filed. Currently \fagin{} performs NO checking on correctness of order. It
is assumed that the input was created by an (as of yet unwritten) utility
that parsed it from a GFF file and carefully validated it.

@<Program input functions@>=
Node * loadNodeList(char * filename){
    FILE * fp = fopen(filename, "rb");
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    if(fp == NULL){
        printf("Cannot open file\n");
        exit(1);
    }
    
    // Define root node
    Node * root = newNode();
    strcpy(root->name, "root");
    root->type = 0;


    // 1. create a flat, doubly-linked list of Node objects
    // the hierarchical level is stored in type {1, 2, ..., nlevels - 1}
    Node * lastnd = root;
    while ((read = getline(&line, &len, fp)) != EOF) {
        Node * nd = newNode();
        sscanf(line, "%*s %d %s", &nd->type, nd->name);
        nd->last = lastnd;
        lastnd->next = nd;
        lastnd = nd;
    }
    lastnd->next = NULL;
    return(root);
}
@

@*2 Loading synteny files.

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



\cwebIndexIntro{%
Below is an index of identifiers. Location of declaration is underlined.
}

@*1 Function prototypes.

@<Function prototypes@>=
void printNode(Node * node);
void recursivePrintNode(Node * node);
Node * newNode();
void freeNode(Node * node);
void freeList(Node * node);
Node * buildTree(Node * root, size_t nlevels);
void wireGenomeTree(Node * node);
Node * loadNodeList(char * filename);
@

\end{document}
