#include "contig.h"

Contig* init_Contig(char* name, size_t size, long length)
{
    Contig* con = (Contig*)malloc(sizeof(Contig));
    con->name   = strdup(name);
    con->length = length;
    con->size   = size;
    con->block  = (Block*)calloc(size, sizeof(Block));
    con->itree  = NULL;
    con->cor[0] = NULL;
    con->cor[1] = NULL;
    con->cor[2] = NULL;
    con->cor[3] = NULL;
    con->cset   = NULL;
    con->ctree  = NULL;
    return con;
}

void free_Contig(Contig* contig)
{
    if (contig != NULL) {
        if (contig->block) {
            free(contig->block);
        }
        if (contig->itree != NULL) {
            free_IntervalTree(contig->itree);
        }
        if (contig->ctree != NULL) {
            free_IntervalTree(contig->ctree);
        }

        while (contig->cset != NULL){
            free_ContiguousSet(contig->cset);
        }

        free(contig->name);
        free(contig);
    }
}

void print_Contig(Contig* contig, bool forward)
{
    fprintf(
        stderr,
        "$ %s size=%lu length=%lu cor=[%zu,%zu,%zu,%zu]\n",
        contig->name,
        contig->size,
        contig->length,
        contig->cor[0]->linkid,
        contig->cor[1]->linkid,
        contig->cor[2]->linkid,
        contig->cor[3]->linkid);
    for (ContiguousSet* cset = contig->cset; cset != NULL; cset = cset->next) {
        fprintf(stderr, "  -- ");
        print_ContiguousSet(contig->cset);
    }
    Block* blk = contig->cor[!forward];
    Corner d = forward ? NEXT_START : NEXT_STOP;
    for (; blk != NULL; blk = blk->cor[d]) {
        print_verbose_Block(blk);
    }
}

ResultContig* init_ResultContig(Contig* contig, IntervalResult* ir, bool is_cset)
{
    ResultContig* rc = (ResultContig*)malloc(sizeof(ResultContig));
    rc->size = ir->iv->size;
    rc->contig = contig;
    rc->inbetween = ir->inbetween;
    rc->leftmost = ir->leftmost;
    rc->rightmost = ir->rightmost;
    if(is_cset){
        rc->cset  = (ContiguousSet**)malloc(rc->size * sizeof(ContiguousSet*));
        rc->block = NULL;
        for (size_t i = 0; i < rc->size; i++) {
            rc->cset[i] = (ContiguousSet*)ir->iv->v[i].link;
        }
    } else {
        rc->block = (Block**)malloc(rc->size * sizeof(Block*));
        rc->cset  = NULL;
        for (size_t i = 0; i < rc->size; i++) {
            rc->block[i] = (Block*)ir->iv->v[i].link;
        }
    }
    return rc;
}

void free_ResultContig(ResultContig* rc)
{
    if (rc != NULL) {
        if (rc->block != NULL) {
            free(rc->block);
        }
        if (rc->cset != NULL) {
            free(rc->cset);
        }
        free(rc);
    }
}

/*
void print_ResultContig(ResultContig * rc)
{
    printf("inbetween=%i leftmost=%i rightmost=%i\n", rc->inbetween, rc->leftmost, rc->rightmost);
    print_Contig(rc->contig, true);
}
*/

void build_block_itree(Contig* con)
{
    // Build itree if necessary
    if (con != NULL && con->itree == NULL){
        Block* blk = con->cor[0];
        // Count the number of blks in the linked list
        // Since blocks can be deleted, we cannot just use con->size
        size_t n = 0;
        for (; blk != NULL; blk = blk->cor[1]) {
            n++;
        }
        IA* ia = init_set_IA(n);
        // Array index for ia->v array of intervals
        size_t i = 0;
        // Reset blk, since was wound to NULL above
        blk = con->cor[0];
        // Map the Block list to the new IA structure
        for (; blk != NULL; blk = blk->cor[1], i++) {
            ia->v[i].start = blk->pos[0];
            ia->v[i].stop  = blk->pos[1];
            ia->v[i].link  = (void*)blk;
        }
        con->itree = build_tree(ia);
    }
}

void build_cset_itree(Contig* con)
{
    // Build ctree if necessary
    if (con != NULL && con->ctree == NULL){
        ContiguousSet * cset = con->cset;
        // rewind cset if necessary
        while(cset->prev != NULL){
            cset = cset->prev;
        }
        size_t n = 0;
        for (; cset != NULL; cset = cset->next) {
            n++;
        }
        // Reset cset, since was wound to NULL above
        cset = con->cset;
        // Map the Block list to the new IA structure
        IA* ia = init_set_IA(n);
        for (size_t i = 0; cset != NULL; cset = cset->next, i++) {
            ia->v[i].start = cset->bounds[0];
            ia->v[i].stop  = cset->bounds[1];
            ia->v[i].link  = (void*)cset;
        }
        con->ctree = build_tree(ia);
    }
}

ResultContig* get_region(Contig * con, long a, long b, bool is_cset)
{

    IntervalTree * tree;

    if(is_cset){
        build_cset_itree(con);
        tree = con->ctree;
    } else {
        build_block_itree(con);
        tree = con->itree;
    }

    // Search itree
    Interval* inv = init_Interval(a, b);
    IntervalResult* res = get_interval_overlaps(inv, tree);
    free(inv);

    // itree returns the flanks for queries that overlap nothing. However, I
    // need all the intervals that overlap these flanks as well.
    IntervalResult* tmp_a = NULL;
    IntervalResult* tmp_b = NULL;

    // I want everything that overlaps the flanks if I am doing a Block search.
    // For ContiguousSet searches, I currently only want overlapping intervals.
    if(!is_cset){
        if (res->inbetween) {
            // If inbetween, itree should have returned the two flanking blocks
            if (res->iv->size == 2) {
                tmp_a = get_interval_overlaps(&res->iv->v[0], con->itree);
                tmp_b = get_interval_overlaps(&res->iv->v[1], con->itree);
                free_IV(res->iv);
                res->iv = tmp_a->iv;
                merge_IV(tmp_a->iv, tmp_b->iv);
                free(tmp_a);
                free_IntervalResult(tmp_b);
            } else {
                fprintf(stderr, "itree is broken, should return exactly 2 intervals for inbetween cases\n");
                exit(EXIT_FAILURE);
            }
        } else if (res->leftmost || res->rightmost) {
            if (res->iv->size == 1) {
                tmp_a = get_interval_overlaps(&res->iv->v[0], con->itree);
                free_IV(res->iv);
                res->iv = tmp_a->iv;
                free(tmp_a);
            } else {
                fprintf(stderr, "itree is broken, should return only 1 interval for left/rightmost cases\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    ResultContig* resultcontig = init_ResultContig(con, res, is_cset);

    free_IntervalResult(res);

    return resultcontig;
}

long count_overlaps(Contig* con, long a, long b)
{
    build_block_itree(con);
    Interval* inv = init_Interval(a, b);
    long count = count_interval_overlaps(inv, con->itree);
    free(inv);
    return count;
}


void merge_doubly_overlapping_blocks(Contig *con)
{
    if (con->cor[0] == NULL ||
        con->cor[1] == NULL ||
        con->cor[2] == NULL ||
        con->cor[3] == NULL)
    {
        fprintf(stderr, "Contig head must be set before merging blocks\n");
        exit(EXIT_FAILURE);
    }

    Block *hi, *lo;

    // iterate up in start-based order
    for(hi = con->cor[0]; hi != NULL; hi = hi->cor[1])
    {
        // look back in stop-based order
        for(lo = hi->cor[3]; lo != NULL; lo = lo->cor[3])
        {
            if(lo->grpid != hi->grpid){
                break;
            }
            if(lo->over->grpid == hi->over->grpid){
                merge_block_a_into_b(lo, hi);
                lo = hi;
            }
        }
    }
}

void sort_blocks(Block** blocks, size_t size, bool by_stop)
{
    if (by_stop) {
        qsort(blocks, size, sizeof(Block*), block_cmp_stop);
    } else {
        qsort(blocks, size, sizeof(Block*), block_cmp_start);
    }
}
