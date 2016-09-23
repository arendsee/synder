#include "block.h"

Block::Block() { }

Block::Block(
    long     t_start,
    long     t_stop,
    double   t_score,
    char     t_strand,
    Feature* t_parent,
    size_t   t_linkid
)
    :
    LinkedInterval(t_parent, t_score, t_strand, t_linkid),
    Interval( t_start, t_stop )
{ }

Block::~Block() { }

void Block::print()
{
    printf("%s\t%zu\t%zu\t%s\t%zu\t%zu\t%c\n",
           parent->name.c_str(),
           pos[0] + Offsets::out_start,
           pos[1] + Offsets::out_stop,
           over->parent->name.c_str(),
           over->pos[0] + Offsets::out_start,
           over->pos[1] + Offsets::out_stop,
           strand
          );
}

void Block::unlink(Block* blk, int u, int d)
{
    if(blk->cor[u] != nullptr) {
        blk->cor[u]->cor[d] = blk->cor[d];
    }
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
    if(a->cor[u] != nullptr) {
        if(a->cor[u] != b) {
            b->cor[u] = a->cor[u];
            a->cor[u]->cor[d] = b;
        }
    } else {
        b->cor[u] = nullptr;
    }
}

void Block::move_b_to_a(Block* a, Block* b, int u, int d)
{

    dissolve_edge(b, u, d);

    replace_edge(a, b, u, d);
    replace_edge(a, b, d, u);

}

void Block::merge_block_a_into_b(Block* a, Block* b)
{
    if(! (a->overlap(b) && a->over->overlap(b->over)) ) {
        fprintf(stderr, "Blocks are not doubly overlapping, I don't know how to merge them\n");
    }

    double olen = a->overlap_length(b);
    double al = a->pos[1] - a->pos[0] + 1;
    double bl = b->pos[1] - b->pos[0] + 1;
    double score =
        ((al - olen) / (al)) * a->score +
        ((bl - olen) / (bl)) * b->score +
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
