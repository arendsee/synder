#include "block.h"

Block::Block() { }

Block::Block(
    long    new_start,
    long    new_stop,
    double  new_score,
    char    new_strand,
    Contig* new_parent,
    size_t  new_linkid
)
{
    parent = new_parent;
    score  = new_score;
    strand = new_strand;
    id     = new_linkid;
    pos[0] = new_start;
    pos[1] = new_stop;
    over   = NULL;
    cor[0] = NULL;
    cor[1] = NULL;
    cor[2] = NULL;
    cor[3] = NULL;
    adj[0] = NULL;
    adj[1] = NULL;
    cnr[0] = NULL;
    cnr[1] = NULL;
    cset   = NULL;
}

Block::~Block() { }

void Block::print()
{
    // printf("%s\t%zu\t%zu\t%s\t%zu\t%zu\t%c",
    //        parent->name.c_str(),
    //        pos[0] + Offsets::out_start,
    //        pos[1] + Offsets::out_stop,
    //        over->parent->name.c_str(),
    //        over->pos[0] + Offsets::out_start,
    //        over->pos[1] + Offsets::out_stop,
    //        strand
    //       );
    // if(cset != NULL) {
    //     printf("\t%zu\n", cset->id);
    // } else {
    //     printf("\t-\n");
    // }
}

void Block::unlink(Block* blk, int u, int d)
{
    // if(blk->cor[u] != NULL) {
    //     blk->cor[u]->cor[d] = blk->cor[d];
    // }
    // if(blk->parent->cor[u] == blk) {
    //     if(blk->cor[d] != NULL) {
    //         blk->parent->cor[u] = blk->cor[d];
    //     } else {
    //         fprintf(stderr, "ERROR: cannot set parent cor to NULL in %s:%d\n", __func__, __LINE__);
    //     }
    // }
}


/** Merge one edge of a into b
 *
 * i can be 0 or 1
 *
 */
void Block::merge_block_a_into_b_edge_(Block* a, Block* b, int i)
{
    assert(i == 0 || i == 1);
    int u = 2 * i +  i;
    int d = 2 * i + !i;
    if(REL_GE(a->pos[i], b->pos[i], i)) {
        move_b_to_a(a, b, u, d);
        b->pos[i] = a->pos[i];
    } else {
        dissolve_edge(a, u, d);
    }
}

void Block::dissolve_edge(Block* blk, int u, int d)
{
    unlink(blk, u, d);
    unlink(blk, d, u);
}

void Block::replace_edge(Block* a, Block* b, int u, int d)
{
    // if(a->cor[u] != NULL) {
    //     if(a->cor[u] != b) {
    //         b->cor[u] = a->cor[u];
    //         a->cor[u]->cor[d] = b;
    //     }
    // } else {
    //     b->cor[u] = NULL;
    // }
    // if(a->parent->cor[u] == a) {
    //     a->parent->cor[u] = b;
    // }
}

void Block::move_b_to_a(Block* a, Block* b, int u, int d)
{

    dissolve_edge(b, u, d);

    replace_edge(a, b, u, d);
    replace_edge(a, b, d, u);

}

void Block::merge_block_a_into_b(Block* a, Block* b)
{
    // if(! (a->overlap(b) && a->over->overlap(b->over)) ) {
    //     fprintf(stderr, "Blocks are not doubly overlapping, I don't know how to merge them\n");
    // }
    //
    // long olen = a->overlap_length(b);
    // long al = a->pos[1] - a->pos[0] + 1;
    // long bl = b->pos[1] - b->pos[0] + 1;
    // double score =
    //     ((double)(al - olen) / (al)) * a->score +
    //     ((double)(bl - olen) / (bl)) * b->score +
    //     olen * (a->score / al + b->score / bl) / 2;
    //
    // b->score       = score;
    // b->over->score = score;
    //
    // merge_block_a_into_b_edge_(a, b, 0);
    // merge_block_a_into_b_edge_(a, b, 1);
    // merge_block_a_into_b_edge_(a->over, b->over, 0);
    // merge_block_a_into_b_edge_(a->over, b->over, 1);
    //
    // memset(a->over, 0, sizeof(Block));
    // memset(a, 0, sizeof(Block));
}
