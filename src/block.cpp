#include "block.h"

/** Initialize values in a block */
void set_Block(
    Block*  block,
    long    start,
    long    stop,
    double   score,
    char    strand,
    Contig* parent,
    Block*  over,
    size_t  linkid
)
{
    block->pos[0] = start;
    block->pos[1] = stop;
    block->score  = score;
    block->strand = strand;
    block->parent = parent;
    block->over   = over;
    block->linkid = linkid;
}

/** Allocate memory for a block and set each field.
 *
 * @param start  query start position
 * @param stop   query stop position
 *
 * @return pointer to new Block
 *
 */
Block* init_Block(long start, long stop)
{
    Block* block  = (Block*)calloc(1, sizeof(Block));
    block->pos[0] = start;
    block->pos[1] = stop;
    block->strand = '.';
    return (block);
}

/** Remove this block from the web of linkage
 *
 * The block is not freed, since it is assumed the memory is in a memory pool
 * managed by a Contig object.
 */
void delete_Block_(Block* block)
{
    int r;
    for(int i = 0; i < 4; i++){
        r = (i % 2 == 0) ? (i + 1) : (i - 1);
        if (block == block->parent->cor[i]){
            block->parent->cor[i] = block->cor[r];
        }
        if (block->cor[i] != NULL){
            block->cor[i]->cor[r] = block->cor[r];
        }
    }

    for(int i = 0; i < 2; i++){
        if (block->cnr[i] != NULL){
            block->cnr[i]->cnr[!i] = block->cnr[!i];
        }
        if (block->adj[i] != NULL) {
            // Normally we link block->adj[1]->adj[0] to block->cor[PREV_STOP], but
            // in the case below, we need to link to NEXT_STOP instead:
            //                    <====>   block->adj[1]
            // i+1       <=====>           block->cor[NEXT_STOP]
            // i         <=====>           block
            // i-1     <====>              block->cor[PREV_STOP]
            r = (block->cor[3*i]->pos[i] < block->adj[i]->pos[!i]) ? 3*i : i+1;
            block->adj[i]->adj[!i] = block->cor[r];
            block->cor[r]->adj[i]  = block->adj[i];
        }
    }

    // wipe memory: unnecessary, but will ease debugging
    memset(block, 0, sizeof(Block));
}

void delete_Block(Block* block)
{
    if(block->over != NULL)
        delete_Block_(block->over);
    if(block != NULL)
        delete_Block_(block);
}

long overlap_length_ll(long a1, long a2, long b1, long b2){
    // If the intervals overlap
    if(a1 <= b2 && b1 <= a2){
        // Find the lower bound of the overlapping region
        long a = a1 > b1 ? a1 : b1; 
        // Find the upper bound of the overlapping region
        long b = a2 > b2 ? b2 : a2;
        // Return the overlapping interval length
        return b - a + 1;
    }
    else{
        return 0;
    }
}

long overlap_length(Block * a, Block * b){
    return overlap_length_ll(
        a->pos[0],
        a->pos[1],
        b->pos[0],
        b->pos[1]
    );
}

void unlink(Block * blk, int u, int d){
    if(blk->cor[u] != NULL){
        blk->cor[u]->cor[d] = blk->cor[d];
    }
    if(blk->parent->cor[u] == blk){
        if(blk->cor[d] != NULL){
            blk->parent->cor[u] = blk->cor[d];
        } else {
            fprintf(stderr, "ERROR: cannot set parent cor to NULL in %s:%d\n", __func__, __LINE__);
        }
    }
}

void dissolve_edge(Block * blk, int u, int d){
    unlink(blk, u, d);
    unlink(blk, d, u);
}

void replace_edge(Block * a, Block * b, int u, int d){
    if(a->cor[u] != NULL){
        if(a->cor[u] != b){
            b->cor[u] = a->cor[u];
            a->cor[u]->cor[d] = b;
        }
    } else {
        b->cor[u] = NULL;
    }
    if(a->parent->cor[u] == a){
        a->parent->cor[u] = b;
    }
}

void move_b_to_a(Block * a, Block * b, int u, int d){

    dissolve_edge(b, u, d);

    replace_edge(a, b, u, d);
    replace_edge(a, b, d, u);

}

/** Merge one edge of a into b
 *
 * i can be 0 or 1
 *
 */
void merge_block_a_into_b_edge_(Block * a, Block * b, int i){
    assert(i == 0 || i == 1);
    int u = 2 * i +  i;
    int d = 2 * i + !i;
    if(REL_GE(a->pos[i], b->pos[i], i)){
        move_b_to_a(a, b, u, d);
        b->pos[i] = a->pos[i];
    } else {
        dissolve_edge(a, u, d);
    }
}

void merge_block_a_into_b(Block * a, Block * b){

    if(! (block_overlap(a, b) && block_overlap(a->over, b->over)) ){
        fprintf(stderr, "Blocks are not doubly overlapping, I don't know how to merge them\n");
    }

    long olen = overlap_length(a, b);
    long al = a->pos[1] - a->pos[0] + 1;
    long bl = b->pos[1] - b->pos[0] + 1;
    double score =
        ((double)(al - olen) / (al)) * a->score +
        ((double)(bl - olen) / (bl)) * b->score +
        olen * (a->score / al + b->score / bl) / 2;

    b->score       = score;
    b->over->score = score;

    merge_block_a_into_b_edge_(a, b, 0);
    merge_block_a_into_b_edge_(a, b, 1);
    merge_block_a_into_b_edge_(a->over, b->over, 0);
    merge_block_a_into_b_edge_(a->over, b->over, 1);

    memset(a->over, 0, sizeof(Block));
    memset(a, 0, sizeof(Block));
}

void print_Block(Block* block)
{
    printf("%s\t%zu\t%zu\t%s\t%zu\t%zu\t%c",
        block->parent->name.c_str(),
        block->pos[0] + Offsets::out_start,
        block->pos[1] + Offsets::out_stop,
        block->over->parent->name.c_str(),
        block->over->pos[0] + Offsets::out_start,
        block->over->pos[1] + Offsets::out_stop,
        block->strand
    );
    if(block->cset != NULL){
        printf("\t%zu\n", block->cset->id);
    } else {
        printf("\t-\n");
    }
}

bool overlap(long a1, long a2, long b1, long b2)
{
    return a1 <= b2 && a2 >= b1;
}

bool block_overlap(Block* a, Block* b)
{
    return a->pos[0] <= b->pos[1] && a->pos[1] >= b->pos[0];
}

/** Compare by Block stop position */
int block_cmp_stop(const void* ap, const void* bp)
{
    Block* a = *(Block**)ap;
    Block* b = *(Block**)bp;
    return (int)(a->pos[1] > b->pos[1]) - (int)(a->pos[1] < b->pos[1]);
}

/** Compare by Block start position */
int block_cmp_start(const void* ap, const void* bp)
{
    Block* a = *(Block**)ap;
    Block* b = *(Block**)bp;
    return (int)(a->pos[0] > b->pos[0]) - (int)(a->pos[0] < b->pos[0]);
}
