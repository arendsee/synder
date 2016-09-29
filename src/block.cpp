#include "block.h"

Block::Block()
{ }

Block::Block(
    long t_start,
    long t_stop,
    double t_score,
    char t_strand,
    Feature* t_parent,
    size_t t_linkid)
    : LinkedInterval(t_parent, t_score, t_strand, t_linkid)
    , Interval(t_start, t_stop)
{ }

Block::~Block() { }

void Block::print()
{
    printf(
        "%s\t%ld\t%ld\t%s\t%ld\t%ld\t%lf\t%c\n",
        parent->name.c_str(),               // query chromosome
        pos[0] + Offsets::out_start,        // query start
        pos[1] + Offsets::out_stop,         // query stop
        over->parent->name.c_str(),         // target chromosome
        over->pos[0] + Offsets::out_start,  // target start
        over->pos[1] + Offsets::out_stop,   // target stop
        score,                              // score
        strand                              // strand
    );
}

void Block::unlink(Block* blk, int u, int d)
{
    if (blk->cor[u] != nullptr) {
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
    int u = 2 * i + i;
    int d = 2 * i + !i;
    if (REL_GE(a->pos[i], b->pos[i], i)) {
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
    if (a->cor[u] != nullptr) {
        if (a->cor[u] != b) {
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
    if (!(a->overlap(b) && a->over->overlap(b->over))) {
        throw "Blocks are not doubly overlapping, I don't know how to merge them";
    }

    double olen = a->overlap_length(b);
    double al = a->pos[1] - a->pos[0] + 1;
    double bl = b->pos[1] - b->pos[0] + 1;
    double score = ((al - olen) / (al)) * a->score +
                   ((bl - olen) / (bl)) * b->score +
                   std::max(a->score / al, b->score / bl) * olen;

#ifdef DEBUG
    std::cerr << "merging a into b: \n"
              << "  -- intermediate values:\n"
              << "  len=(a:" << al << ",b:" << bl << ") overlap_len = " << olen << " scores=(a:" << a->score << ",b:" << b->score << ")\n"
              << "  -- original block\n"
              << "  Q:a = (score:" << a->score       << " id:" << a->id       << " bounds=[" << a->pos[0]       << "," << a->pos[1]       << ")\n"
              << "  Q:b = (score:" << b->score       << " id:" << b->id       << " bounds=[" << b->pos[0]       << "," << b->pos[1]       << ")\n"
              << "  T:a = (score:" << a->over->score << " id:" << a->over->id << " bounds=[" << a->over->pos[0] << "," << a->over->pos[1] << ")\n"
              << "  T:b = (score:" << b->over->score << " id:" << b->over->id << " bounds=[" << b->over->pos[0] << "," << b->over->pos[1] << ")\n";
#endif

    b->score = score;
    b->over->score = score;

    merge_block_a_into_b_edge_(a, b, 0);
    merge_block_a_into_b_edge_(a, b, 1);
    merge_block_a_into_b_edge_(a->over, b->over, 0);
    merge_block_a_into_b_edge_(a->over, b->over, 1);

#ifdef DEBUG
    std::cerr << "  -- merged blocks\n"
              << "  Q:ab = (score:" << b->score       << " id:" << b->id       << " bounds=[" << b->pos[0]       << "," << b->pos[1] << ")\n"
              << "  T:ab = (score:" << b->over->score << " id:" << b->over->id << " bounds=[" << b->over->pos[0] << "," << b->over->pos[1] << ")\n";
#endif

    // Declare these blocks broken. A null `over` tags these blocks for
    // exclusion and will break asserts in Genome::validate
    a->over->over = nullptr;
    a->over = nullptr;
}
