#include "block.h"

/** Initialize values in a block */
void set_Block(
    Block*  block,
    long    start,
    long    stop,
    float   score,
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
    Block* block = (Block*)calloc(1, sizeof(Block));
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


void merge_block_a_into_b(Block * a, Block * b){

    // 1. reasign b->score to weighted average of a->score and b->score

    long olen = overlap_length(a, b);
    long al = a->pos[1] - a->pos[0] + 1;
    long bl = b->pos[1] - b->pos[0] + 1;
    float score =
        ((float)(al - olen) / (al)) * a->score +
        ((float)(bl - olen) / (bl)) * b->score +
        olen * (a->score / al + b->score / bl) / 2;

    b->score = score;
    b->over->score = score;
        
    // 2. reasign b->pos to union of a->pos and b->pos

    b->pos[0] = b->pos[0] < a->pos[0] ? b->pos[0] : a->pos[0];
    b->pos[1] = b->pos[1] > a->pos[1] ? b->pos[1] : a->pos[1];
    b->over->pos[0] = b->over->pos[0] < a->over->pos[0] ? b->over->pos[0] : a->over->pos[0];
    b->over->pos[1] = b->over->pos[1] > a->over->pos[1] ? b->over->pos[1] : a->over->pos[1];

    // 3. delete a, this will automatically relink

    delete_Block(a);
}

void print_Block(Block* block)
{
    printf("%s\t%zu\t%zu\t%s\t%zu\t%zu\t%c",
        block->parent->name,
        block->pos[0] + global_out_start,
        block->pos[1] + global_out_stop,
        block->over->parent->name,
        block->over->pos[0] + global_out_start,
        block->over->pos[1] + global_out_stop,
        block->strand
    );
    if(block->cset != NULL){
        printf("\t%zu\n", block->cset->id);
    } else {
        printf("\t-\n");
    }
}

// if block is NULL, "-"
// else string of size_t
void link_str(char * s, Block * b){
    if(b == NULL) {
        strcpy(s, "-");
    }
    else{
        sprintf(s, "%zu", b->linkid);
    }
}
void print_verbose_Block_(Block* block, char side)
{
    fprintf(
        stderr,
        "  %c-%zu parent=%s pos(%zu, %zu, %c) grpid=%zu addr=%p\n",
        side,
        block->linkid,
        block->parent->name,
        block->pos[0],
        block->pos[1],
        block->strand,
        block->grpid,
        block
    );
    char c0[32], c1[32], c2[32], c3[32], a0[32], a1[32], r0[32], r1[32];
    link_str(c0, block->cor[0]);
    link_str(c1, block->cor[1]);
    link_str(c2, block->cor[2]);
    link_str(c3, block->cor[3]);
    link_str(a0, block->adj[0]);
    link_str(a1, block->adj[1]);
    link_str(r0, block->cnr[0]);
    link_str(r1, block->cnr[1]);
    fprintf(
        stderr,
        "  %clinks\n    * cor=[%s,%s,%s,%s]\n    * adj=[%s,%s]\n    * cnr=[%s,%s]\n",
        side, c0, c1, c2, c3, a0, a1, r0, r1
    );
}
void print_verbose_Block(Block* block)
{
    fprintf(
        stderr,
        "$ setid=%zu score=%f\n",
        block->cset->id,
        block->score
    );
    print_verbose_Block_(block, 'Q');
    print_verbose_Block_(block->over, 'T');
}

/**
 * Determine whether interval (a,b) overlaps interval (c,d)
 *
 * @param a1 start of first interval
 * @param a2 stop of first interval
 * @param b1 start of second interval
 * @param b2 stop of second interval
 *
 * @return TRUE if the intervals overlap
 */
bool overlap(long a1, long a2, long b1, long b2)
{
    return a1 <= b2 && a2 >= b1;
}

/**
 * Determine whether two Blocks overlap 
 *
 * @return TRUE if they overlap 
 */
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

long get_set_bound(Block* blk, Direction d)
{
    if (blk->cnr[d] != NULL) {
        return get_set_bound(blk->cnr[d], d);
    }
    return blk->pos[d];
}
