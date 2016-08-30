fs next

# break -function load_Synmap
# break -function link_block_corners
# break -function set_contig_corners
# break -function set_overlap_group
# break -function link_adjacent_blocks

break -function link_contiguous_blocks
break -function find_search_intervals

# set print repeats 100
# 
# define ps
#     dont-repeat
#     print "--------------------------------------------"
#     if $argc == 0
#         print *syn
#         print syn->genome[0]
#         print syn->genome[1]
#     end
#     if $argc == 1
#         print *syn->genome[$arg0]
#     end
#     if $argc == 2
#         print *syn->genome[$arg0]->contig[$arg1]
#     end
# end
# document ps
#     Print synmap, args = [0,1,2]
# end
# 
# define pb
#     dont-repeat
#     print "--------------------------------------------"
#     print *blk
# end
# document pb
#     Print block, no args
# end
# 
# define pbp
#     dont-repeat
#     print "--------------------------------------------"
#     print *blk->parent
# end
# document pbp
#     Print block parent
# end
# 
# define pbo
#     dont-repeat
#     print "--------------------------------------------"
#     print *blk->over
# end
# document pbo
#     Print block's homolog
# end
# 
# define pb-cor
#     dont-repeat
#     print "--------------------------------------------"
#     print "prev start"
#     print *blk->cor[0]
#     print "next start"
#     print *blk->cor[1]
#     print "prev stop"
#     print *blk->cor[2]
#     print "next stop"
#     print *blk->cor[3]
# end
# document pb-cor
#     Print block's 4 corners
# end
# 
# define pb-adj
#     dont-repeat
#     print "--------------------------------------------"
#     print *blk->adj[0]
#     print *blk->adj[1]
# end
# document pb-adj
#     Print block's adj pair
# end
# 
# define pb-cnr
#     dont-repeat
#     print "--------------------------------------------"
#     print *blk->cnr[0]
#     print *blk->cnr[1]
# end
# document pb-adj
#     Print block's contiguous neighbors
# end
# 
# define psv
#     dont-repeat
#     print "--------------------------------------------"
#     if $argc == 1
#         if $arg0 == 0
#             print *syn
#         end
#         if $arg0 == 1
#             print *syn
#             print *syn->genome[0]
#             print *syn->genome[1]
#         end
# 
#         if $arg0 == 2
#             print *syn->genome[0]->contig[0]
#             print *syn->genome[1]->contig[0]
#         end
# 
#     end
# end
# document psv
#     Print synmap recursively args = [012]
# end
# 
# define pcb
#    dont-repeat
#    print "--------------------------------------------"
#    print "---print all blocks---"
#    set $i = (int)syn->genome[$arg0]->size - 1
#    while $i >= 0
#        set $j = (int)syn->genome[$arg0]->contig[$i]->size - 1
#        while $j >= 0
#            print syn->genome[$arg0]->contig[$i]->block[$j]
#            set $j= $j - 1
#        end
#        set $i= $i - 1
#    end
# end
# document pcb
#     Print all blocks in genome, args=[01]
# end
# 
# define pc
#     dont-repeat
#     print "--------------------------------------------"
#     print *con
#     print *con->cor[0]
#     print *con->cor[1]
# end
# document pc
#     Print contig corners 0 and 1
# end
# 
# define pnode
#     dont-repeat
#     print "--------------------------------------------"
#     if $argc == 0
#         print *node 
#         print  node->blk->parent->name
#         print *node->blk
#         print "over"
#         print  node->blk->over->parent->name
#         print *node->blk->over
#     end
#     if $argc == 1
#         if $arg0 == 1
#             print *node->down 
#             print  node->down->blk->parent->name
#             print *node->down->blk
#             print "over"
#             print  node->blk->over->parent->name
#             print *node->blk->over
#         end
#         if $arg0 == 2
#             print *node->down->down 
#             print  node->down->down->blk->parent->name
#             print *node->down->down->blk
#             print "over"
#             print  node->down->down->blk->over->parent->name
#             print *node->down->down->blk->over
#         end
#         if $arg0 == 3
#             print *node->down->down->down 
#             print  node->down->down->down->blk->parent->name
#             print *node->down->down->down->blk
#             print "over"
#             print  node->down->down->down->blk->over->parent->name
#             print *node->down->down->down->blk->over
#         end
#     end
# end
# document pnode
#     Recursively print Node to [0123] levels down
# end
