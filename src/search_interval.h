#ifndef __SEARCH_INTERVAL_H__
#define __SEARCH_INTERVAL_H__

class SearchInterval : Interval<SearchInterval>
{
private:
    double score;
    char* seqname;

    // TODO create SI class
    typedef struct SI_Bound {
        long bound;
        int flag;
    } SI_Bound;

    SI_Bound* get_si_bound(
        long q,
        Block* blk_bounds[2],
        Direction d,
        bool inverted
    );

    // TODO remove the local LL
    // A local utility structure used filter and store contiguous sets
    typedef struct CSList {
        struct CSList* next;
        Block* bound[2];
        ContiguousSet* cset;
    } CSList;
    CSList* init_empty_CSList();
    CSList* init_CSList(Block* blk);
    void add_blk_CSList(CSList* cslist, Block* blk);
    void add_cset_CSList(CSList* cslist, ContiguousSet* cset, long bounds[2]);
    void free_CSList(CSList* cslist);

    double calculate_score(Bound& bound, Block* blk)
    double flank_area(long near, long far, double k);

    int get_flag(SI_Bound* br[2]);

public:
    SearchInterval(Bound bound, char* seqname);

    void print();

}



#endif
