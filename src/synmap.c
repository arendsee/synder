#include "global.h"
#include "synmap.h"
#include "contiguous_set.h"


Synmap* init_Synmap()
{
    Synmap* syn = (Synmap*)malloc(sizeof(Synmap));
    syn->size = 2;
    syn->genome = (Genome**)malloc(2 * sizeof(Genome*));
    return (syn);
}

void free_Synmap(Synmap* synmap)
{
    if (synmap != NULL) {
        free_Genome(SG(synmap, 0));
        free_Genome(SG(synmap, 1));
        free(synmap->genome);
        free(synmap);
    }
}

void print_Synmap(Synmap* synmap, bool forward)
{
    // only print the query Genome, the print_verbose_Block function will print
    // the target information as well
    fprintf(
        stderr,
        "--- Query=(%s, %zu), Target=(%s, %zu)\n",
        SG(synmap, 0)->name,
        SG(synmap, 0)->size,
        SG(synmap, 1)->name,
        SG(synmap, 1)->size
    );
    fprintf(stderr, "---------------------------------------------------------\n");
    print_Genome(SG(synmap, 0), forward);
    // print_Genome(SG(synmap, 1), forward);
}

void link_block_corners(Synmap* syn)
{
    Contig* con;
    size_t N;
    for (size_t g = 0; g <= 1; g++) {
        for (size_t c = 0; c < SG(syn, g)->size; c++) {
            con = SGC(syn, g, c);
            N = con->size;

            Block** blocks = (Block**)malloc(N * sizeof(Block*));
            for (size_t i = 0; i < N; i++) {
                blocks[i] = &con->block[i];
            }

            // sort by stop
            sort_blocks(blocks, N, true);
            for (size_t i = 0; i < N; i++) {
                blocks[i]->cor[PREV_STOP] = (i == 0) ? NULL : blocks[i - 1];
                blocks[i]->cor[NEXT_STOP] = (i == N - 1) ? NULL : blocks[i + 1];
            }
            // sort by start
            sort_blocks(blocks, N, false);
            for (size_t i = 0; i < N; i++) {
                blocks[i]->cor[PREV_START] = (i == 0) ? NULL : blocks[i - 1];
                blocks[i]->cor[NEXT_START] = (i == N - 1) ? NULL : blocks[i + 1];
            }

            free(blocks);
        }
    }
}

void set_contig_corners(Synmap* syn)
{
    Contig* con;
    for (size_t g = 0; g < 2; g++) {
        for (size_t c = 0; c < SG(syn, g)->size; c++) {
            con = SGC(syn, g, c);
            con->cor[1] = &con->block[0];
            con->cor[0] = &con->block[0];
            con->cor[2] = &con->block[con->size - 1];
            con->cor[3] = &con->block[con->size - 1];
            while (con->cor[1]->cor[PREV_STOP] != NULL) {
                con->cor[1] = con->cor[1]->cor[PREV_STOP];
            }
            while (con->cor[0]->cor[PREV_START] != NULL) {
                con->cor[0] = con->cor[0]->cor[PREV_START];
            }
            while (con->cor[2]->cor[NEXT_START] != NULL) {
                con->cor[2] = con->cor[2]->cor[NEXT_STOP];
            }
            while (con->cor[3]->cor[NEXT_STOP] != NULL) {
                con->cor[3] = con->cor[3]->cor[NEXT_STOP];
            }
        }
    }
}

void set_overlap_group(Synmap* syn)
{

    // Holds current overlapping group id
    size_t grpid = 1;
    // Needed for determining overlaps and thus setids
    long maximum_stop = 0;
    // The stop position of the current interval
    long this_stop = 0;
    // Current Block in linked list
    Block* blk;

    // Loop through target and query genomes
    // g := genome id (0 is query, 1 is target)
    for (size_t g = 0; g <= 1; g++) {
        // Loop through each contig in the query genome
        // i := contig id
        for (size_t i = 0; i < SG(syn, g)->size; i++) {
            maximum_stop = 0;
            // Loop through each Block in the linked list
            blk = SGC(syn, g, i)->cor[0];
            for (; blk != NULL; blk = blk->cor[1]) {
                this_stop = blk->pos[1];
                // If the start is greater than the maximum stop, then the block is in
                // a new adjacency group. For this to work, Contig->block must be
                // sorted by start. This sort is performed in build_tree.
                if (blk->pos[0] > maximum_stop) {
                    grpid++;
                }
                if (this_stop > maximum_stop) {
                    maximum_stop = this_stop;
                }
                blk->grpid = grpid;
            }
            // increment to break adjacency between contigs and genomes
            grpid++;
        }
    }
}

void link_adjacent_blocks_directed(Contig* con, Direction d)
{
    // In diagrams:
    // <--- indicates a hi block
    // ---> indicates a lo block
    // All diagrams and comments relative to the d==HI direction

    if (con->cor[d] == NULL || con->cor[!d] == NULL) {
        fprintf(stderr, "Contig head must be set before link_adjacent_blocks is called\n");
        exit(EXIT_FAILURE);
    }

    // Transformed indices for Block->cor
    // ----------------------------------
    // a   b          c   d
    // <--->   <--->  <--->
    //   a - previous element by start
    //   b - previous element by stop
    //   c - next element by start
    //   d - next element by stop
    int idx_a = (!d * 2) + !d; // - 0
    //int idx_b = (d  * 2) + !d; // - 2
    int idx_c = (!d * 2) + d;  // - 1
    int idx_d = (d  * 2) + d;  // - 3

    Block *lo, *hi;

    lo = con->cor[idx_c]; // first element by stop
    hi = con->cor[idx_a]; // first element by start

    while (hi != NULL) {

        //       --->
        // <---
        // OR
        // --->
        //   <---
        // This should should occur only at the beginning
        if (REL_LE(hi->pos[!d], lo->pos[d], d)) {
            hi->adj[!d] = NULL;
            hi = hi->cor[idx_c];
        }

        //  lo     next
        // ---->  ---->
        //               <---
        // If next is closer, and not overlapping the hi, increment lo
        // You increment lo until it is adjacent to the current hi
        else if (REL_LT(lo->cor[idx_d]->pos[d], hi->pos[!d], d)) {
            lo = lo->cor[idx_d];
        }

        // --->
        //      <---
        // The current lo is next to, and not overlapping, current hi
        else {
            hi->adj[!d] = lo;
            hi = hi->cor[idx_c];
        }
    }
}
void link_adjacent_blocks(Synmap* syn)
{
    for (size_t genid = 0; genid <= 1; genid++) {
        for (size_t conid = 0; conid < SG(syn, genid)->size; conid++) {
            link_adjacent_blocks_directed(SGC(syn, genid, conid), HI);
            link_adjacent_blocks_directed(SGC(syn, genid, conid), LO);
        }
    }
}

// ---- A local utility structure used to build contiguous sets ----
typedef struct Node {
    struct Node * down;
    ContiguousSet * cset;
} Node;
Node* init_node(Block* blk)
{
    Node* node = (Node*)malloc(sizeof(Node));
    node->down = NULL;
    node->cset = init_ContiguousSet(blk);
    return (node);
}
void remove_node(Node* node)
{
    if (node->down != NULL) {
        Node* tmp = node->down;
        memcpy(node, node->down, sizeof(Node));
        free(tmp);
    } else {
        node->cset = NULL;
    }
}
void free_node(Node* node)
{
    if (node->down != NULL)
        free_node(node->down);
    free(node);
}

void link_contiguous_blocks(Synmap* syn, long k)
{

    Block *blk;
    Node* node;
    Node* root;
    bool block_added;

    for (size_t i = 0; i < SG(syn, 0)->size; i++) {

        // Initialize the first block in the scaffold
        blk  = SGC(syn, 0, i)->cor[0];
        node = init_node(blk);
        root = node;
        for (blk = blk->cor[1]; blk != NULL; blk = blk->cor[1]) {
            while (true) {

                block_added = add_block_to_ContiguousSet(node->cset, blk, k);

                // if block has joined a set
                if (block_added){
                    // break and process next block
                    break;
                }
                // if we reach the bottom of the Node list
                else if (node->down == NULL) {
                    // begin new contiguous set
                    node->down = init_node(blk);
                    break;
                }
                // if this node is done
                else if (strictly_forbidden(node->cset->ends[1], blk, k)) {
                    // remove it, exposing the next Node 
                    remove_node(node);
                }
                // otherwise
                else {
                    // descend to the next Node
                    node = node->down;
                }
            }
            node = root;
        }
        free_node(root);
    }

    // TODO: remove this, I don't really need the id
    ContiguousSet * cset;
    size_t setid = 0;
    for (size_t i = 0; i < SG(syn, 0)->size; i++) {
        cset = SGC(syn, 0, i)->cset;
        // rewind - will remove it in production code
        while(cset->prev != NULL){
            cset = cset->prev;
        }
        SGC(syn, 0, i)->cset = cset;
        for(; cset != NULL; cset = cset->next){
            setid++;
            cset->id       = setid;
            cset->over->id = setid;
        }
    }

}

void validate_synmap(Synmap* syn)
{
    #define ASSERT(t) if(!(t)){                                            \
                        is_good=false;                                     \
                        if(blk != NULL)                                    \
                            fprintf(stderr, "%zu:%zu ", gid, blk->linkid); \
                        fprintf(stderr, "Assert failed: `" #t "`\n");      \
                      }

    size_t  gid     = 0;
    size_t  cid     = 0;
    size_t  nblks   = 0;
    Contig* con     = NULL;
    Block*  blk     = NULL;
    bool    is_good = true;

    ASSERT(syn->size == 2)
    for (gid = 0; gid < syn->size; gid++) {
        for (cid = 0; cid < SG(syn, gid)->size; cid++) {
            con = SGC(syn, gid, cid);
            blk = con->cor[0];
            nblks = 0;
            for (; blk != NULL; blk = blk->cor[1]) {
                if (!(blk->pos[1] < con->length)) {
                    fprintf(stderr,
                        "WARNING: stop greater than contig length: %zu vs %zu\n",
                        blk->pos[1], con->length);
                }

                nblks++;

                ASSERT(blk->over != NULL);
                ASSERT(blk->cset != NULL);
                ASSERT(blk == blk->over->over);
                ASSERT(blk->cset->id == blk->over->cset->id);
                ASSERT(blk->cset->over == blk->over->cset);
                ASSERT(blk->score == blk->over->score);

                // grpid == 0 only if unset
                ASSERT(blk->grpid != 0);

                if (blk->cnr[1] != NULL) {
                    ASSERT(blk->grpid != blk->cnr[1]->grpid);
                    ASSERT(blk->cset == blk->cnr[1]->cset);
                    ASSERT(blk->cnr[1]->over->cnr[0] != NULL);
                    ASSERT(blk->cnr[1]->over->cnr[0]->over == blk);
                }
                if (blk->cor[NEXT_START] != NULL) {
                    ASSERT(blk->pos[0] <= blk->cor[NEXT_START]->pos[0]);
                    ASSERT(blk->cor[NEXT_START]->cor[PREV_START] != NULL);
                    ASSERT(blk == blk->cor[NEXT_START]->cor[PREV_START]);
                }
                if (blk->cor[NEXT_STOP] != NULL) {
                    ASSERT(blk->pos[1] <= blk->cor[NEXT_STOP]->pos[1]);
                    ASSERT(blk->cor[NEXT_STOP]->cor[PREV_STOP] != NULL);
                    ASSERT(blk == blk->cor[NEXT_STOP]->cor[PREV_STOP]);
                }
                if (blk->cor[PREV_START] != NULL) {
                    ASSERT(blk->cor[PREV_START]->cor[NEXT_START] != NULL);
                }
                if (blk->cor[PREV_STOP] != NULL) {
                    ASSERT(blk->cor[PREV_STOP]->cor[NEXT_STOP] != NULL);
                }
            }
            // blocks may be deleted, so nblks == con->size may not hold
            // however no blocks should ever be added
            ASSERT(nblks <= con->size);
        }
    }
    if(! is_good){
        print_Synmap(syn, true);
        exit(EXIT_FAILURE);
    }
    #undef vprint
}
