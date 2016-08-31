#include "synmap.h"

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
    fprintf(stderr, "Target contigs:\n");
    print_Genome(SG(synmap, 1), forward, false);
    fprintf(stderr, "---------------------------------------------------------\n");
    fprintf(stderr, "Query contigs and blocks:\n");
    print_Genome(SG(synmap, 0), forward, true);
}

void dump_blocks(Synmap* synmap)
{
    for(size_t i = 0; i < SG(synmap, 0)->size; i++){
        Block * blk = SGC(synmap, 0, i)->cor[0];
        for(; blk != NULL; blk = blk->cor[1]){ 
            print_Block(blk);   
        } 
    }    
}

void dump_verbose_block_mem(Synmap* synmap)
{
    for(size_t i = 0; i < SG(synmap, 0)->size; i++){
        if(i == 1){
            fprintf(stderr, "\n");
        }
        for(int j = 0; j < SGC(synmap, 0, i)->size; j++){
            Block * blk = &SGC(synmap, 0, i)->block[j];
            if(blk != NULL){
                print_verbose_Block(blk);   
            } else {
                fprintf(stderr, "NULL\n");
            }
        }
    }    
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
                blocks[i]->cor[2] = (i == 0)     ? NULL : blocks[i - 1];
                blocks[i]->cor[3] = (i == N - 1) ? NULL : blocks[i + 1];
            }
            // sort by start
            sort_blocks(blocks, N, false);
            for (size_t i = 0; i < N; i++) {
                blocks[i]->cor[0] = (i == 0)     ? NULL : blocks[i - 1];
                blocks[i]->cor[1] = (i == N - 1) ? NULL : blocks[i + 1];
            }

            free(blocks);
        }
    }
}

void set_contig_corners(Synmap* syn)
{
    Contig* con;
    int k;
    for (size_t g = 0; g < 2; g++) {
        for (size_t c = 0; c < SG(syn, g)->size; c++) {
            con = SGC(syn, g, c);
            for(int i = 0; i < 4; i++){
                k = i % 2 == 0 ? 0 : con->size - 1;
                con->cor[i] = &con->block[k];
                while (con->cor[i]->cor[i] != NULL) {
                    con->cor[i] = con->cor[i]->cor[i];
                }
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

    if (con->cor[0] == NULL || con->cor[1] == NULL || con->cor[2] == NULL || con->cor[3] == NULL) {
        fprintf(stderr, "Contig head must be set before link_adjacent_blocks is called\n");
        fprintf(stderr, "genome=(%s) contig=(%s)\n", con->parent->name, con->name);
        exit(EXIT_FAILURE);
    }

    // Transformed indices for Block->cor and Contig->cor
    int idx_a = (!d * 2) + !d; // - 0 previous/first element by start
    int idx_b = (!d * 2) +  d; // - 1 next/last element by start
    int idx_c = ( d * 2) + !d; // - 2 previous/first element by stop
    int idx_d = ( d * 2) +  d; // - 3 next/last element by stop

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
            hi = hi->cor[idx_b];
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
            hi = hi->cor[idx_b];
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

void merge_all_doubly_overlapping_blocks(Synmap * syn)
{
    Contig * con;
    for(size_t i = 0; i < SG(syn, 0)->size; i++){
        con = SGC(syn, 0, i);
        merge_doubly_overlapping_blocks(con);
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

    Contig * con;
    Block *blk;
    Node* node;
    Node* root;
    bool block_added;

    for (size_t i = 0; i < SG(syn, 0)->size; i++) {
        // The current Contig
        con = SGC(syn, 0, i);
        // Initialize the first block in the scaffold
        blk  = con->cor[0];
        // A container for the nascent ContiguousSets
        node = init_node(blk);
        // For returning after descent, since Node has no prev
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

    size_t setid = 0;
    ContiguousSet * cset;
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < SG(syn, i)->size; j++){
           con = SGC(syn, i, j); 
            // rewind - TODO - build the csets so this isn't necessary
            while(con->cset->prev != NULL){
                con->cset = con->cset->prev;
            }
            if(i == 0){
                // TODO - remove setids. I don't use them for anything but debugging
                for(cset = con->cset; cset != NULL; cset = cset->next){
                    setid++;
                    cset->id = setid;
                    cset->over->id = setid;
                }
            }
        }
    }
}

void validate_synmap(Synmap* syn)
{
    #define ASSERT_CON(t, con)                               \
        if(!(t)){                                            \
          is_good=false;                                     \
          fprintf(stderr, "Assert failed: `" #t "`\n");      \
          if(blk != NULL)                                    \
            print_Contig(con, false, false);                 \
        }

    #define ASSERT_BLK(t, blk)                               \
        if(!(t)){                                            \
          is_good=false;                                     \
          fprintf(stderr, "Assert failed: `" #t "`\n");      \
          if(blk != NULL){print_verbose_Block(blk);}         \
          else {fprintf(stderr, "NULL\n");}                  \
        }

    size_t  gid     = 0;
    size_t  cid     = 0;
    size_t  nblks   = 0;
    Contig* con     = NULL;
    Block*  blk     = NULL;
    bool    is_good = true;

    assert(syn->size == 2);
    for (gid = 0; gid < syn->size; gid++) {
        for (cid = 0; cid < SG(syn, gid)->size; cid++) {
            con = SGC(syn, gid, cid);

            ASSERT_CON(con->cor[0] != NULL, con);
            ASSERT_CON(con->cor[1] != NULL, con);
            ASSERT_CON(con->cor[2] != NULL, con);
            ASSERT_CON(con->cor[3] != NULL, con);

            blk = con->cor[0];
            nblks = 0;
            for (; blk != NULL; blk = blk->cor[1]) {
                if (!(blk->pos[1] < con->length)) {
                    fprintf(stderr,
                        "WARNING: stop greater than contig length: %zu vs %zu\n",
                        blk->pos[1], con->length);
                }

                nblks++;

                ASSERT_BLK(blk->over != NULL, blk);
                ASSERT_BLK(blk->cset != NULL, blk);
                ASSERT_BLK(blk == blk->over->over, blk);
                ASSERT_BLK(blk->cset->id == blk->over->cset->id, blk);
                ASSERT_BLK(blk->cset->over == blk->over->cset, blk);
                ASSERT_BLK(blk->score == blk->over->score, blk);

                ASSERT_BLK(blk->pos[0] >= con->cor[0]->pos[0], blk);
                ASSERT_BLK(blk->pos[0] <= con->cor[1]->pos[0], blk);
                ASSERT_BLK(blk->pos[1] >= con->cor[2]->pos[1], blk);
                ASSERT_BLK(blk->pos[1] <= con->cor[3]->pos[1], blk);

                if (blk->cor[0] != NULL){
                    ASSERT_BLK(blk->pos[0] >= blk->cor[0]->pos[0], blk);
                    ASSERT_BLK(blk->cor[0]->cor[1]->pos[0] == blk->pos[0], blk);
                }
                if (blk->cor[1] != NULL){
                    ASSERT_BLK(blk->pos[0] <= blk->cor[1]->pos[0], blk);
                    ASSERT_BLK(blk->cor[1]->cor[0]->pos[0] == blk->pos[0], blk);
                }
                if (blk->cor[2] != NULL){
                    ASSERT_BLK(blk->pos[1] >= blk->cor[2]->pos[1], blk);
                    ASSERT_BLK(blk->cor[2]->cor[3]->pos[1] == blk->pos[1], blk);
                }
                if (blk->cor[3] != NULL){
                    ASSERT_BLK(blk->pos[1] <= blk->cor[3]->pos[1], blk);
                    ASSERT_BLK(blk->cor[3]->cor[2]->pos[1] == blk->pos[1], blk);
                }

                // grpid == 0 only if unset
                ASSERT_BLK(blk->grpid != 0, blk);

                for(int i = 0; i < 2; i++){
                    if (blk->cnr[i] != NULL) {
                        ASSERT_BLK(blk->grpid != blk->cnr[i]->grpid, blk);
                        ASSERT_BLK(blk->cset == blk->cnr[i]->cset, blk);
                        ASSERT_BLK(blk->cnr[i]->over->cnr[!i]->over == blk, blk);
                    }
                }

            }
            // blocks may be deleted, so nblks == con->size may not hold
            // however no blocks should ever be added
            ASSERT_CON(nblks <= con->size, con);
        }
    }
    if(! is_good){
        print_Synmap(syn, true);
        exit(EXIT_FAILURE);
    }
    #undef ASSERT_BLK
    #undef ASSERT_CON
}
