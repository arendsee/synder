#include "many_contiguous_sets.h"

void ManyContiguousSets::link_contiguous_blocks(ManyBlocks& many_blks, long k)
{
    auto iter = inv.rbegin();
    Block* blk = many_blks.front()
    for (; blk != NULL; blk = blk->next()) {
        for (iter = inv.rbegin(); ; iter++) {
            if (iter == inv.end()) {
                // if block fits in no set, create a new one
                inv.push_back(new ContiguousSet(blk));
                break;
            }
            // if block successfully joins a set
            else if ((*iter)->add_block(blk, k)) {
                break;
            }
            // if set terminates
            else if (strictly_forbidden((*iter)->ends[1], blk, k)) {
                inv.push_back(new ContiguousSet(blk));
                break;
            }
        }
    }

    size_t setid = 0;
    for(auto c : inv){
        setid++; 
        c->id = setid;
        c->over->id = setid;
    }

}
