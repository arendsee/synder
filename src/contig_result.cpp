#include "contig_result.h"

ResultContig* init_ResultContig(Contig* contig, IntervalResult* ir, bool is_cset)
{
    ResultContig* rc = (ResultContig*)malloc(sizeof(ResultContig));
    rc->size = ir->iv.size();
    rc->contig = contig;
    rc->inbetween = ir->inbetween;
    rc->leftmost = ir->leftmost;
    rc->rightmost = ir->rightmost;
    if (is_cset)
    {
        rc->cset  = (ContiguousSet**)malloc(rc->size * sizeof(ContiguousSet*));
        rc->block = NULL;
        for (size_t i = 0; i < rc->size; i++)
        {
            rc->cset[i] = (ContiguousSet*)ir->iv[i]->link;
        }
    }
    else
    {
        rc->block = (Block**)malloc(rc->size * sizeof(Block*));
        rc->cset  = NULL;
        for (size_t i = 0; i < rc->size; i++)
        {
            rc->block[i] = (Block*)ir->iv[i]->link;
        }
    }
    return rc;
}

void free_ResultContig(ResultContig* rc)
{
    if (rc != NULL)
    {
        if (rc->block != NULL)
        {
            free(rc->block);
        }
        if (rc->cset != NULL)
        {
            free(rc->cset);
        }
        free(rc);
    }
}
