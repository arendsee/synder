#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "contiguous.h"
#include "contig.h"
#include "synmap.h"


ContiguousMap * init_contiguous_map(size_t size){
	ContiguousMap* cmap = (ContiguousMap*)malloc(sizeof(ContiguousMap));
    cmap->size = size;
    cmap->map = (ContiguousNode**)malloc(size * sizeof(ContiguousNode*));
	
	return cmap;
}

void contiguous_query(Synmap * syn, FILE * intfile){
	// Ensure blocks are sorted
	sort_all_contigs(syn);
	
	// count total number of unique block-block pairs for hashmap
	ContiguousMap *cmap= populate_contiguous_map(syn);
    
	char seqname[128];
    int chrid, start, stop,flag;
    Contig * contigs;
    Contig * tcon;
    Block * qblk;
	Block * endblk;
	ContiguousNode * qnode;
    Block * tblk;
	Block twoblk;
    bool missing;
    while ((fscanf(intfile,
                   "%d %*s %*s %d %d %*c %*c %*s %s\n",
                   &chrid, &start, &stop, seqname)) != EOF)
    {
        contigs = get_region(SGC(syn, 0, chrid), start, stop);
		uint32_t region[2]={0,0};
        missing = false;

	    if(contigs->size == 2){
 	      twoblk.start = start;
	        twoblk.stop = stop;
    	    if(!((CB(contigs, 0) && block_overlap(CB(contigs, 0), &twoblk)) ||
        	     (CB(contigs, 1) && block_overlap(CB(contigs, 1), &twoblk))
            	)){
            		missing = true;
        	}
     	}

        for(int i = 0; i < contigs->size; i++){
            qblk = contigs->block[i];
            if(qblk){
                tblk = QT_SGCB(syn, qblk);
				endblk = cmap->map[tblk->linkid]->match;
				if(endblk->oblkid < region[0]){region[0]=endblk->oblkid;}
				if(endblk->oblkid > region[1]){region[1]=endblk->oblkid;}
        	}
		}
		for(int i= region[0]; i<= region[1];i++){
		   /*
 			* Flag -> keeps track of edge reliability 
 			* 0 = determinable boundries
 			* 1 = Left hand  undeterminable (case F)
 			* 2 = Right hand undertrminable (caseF)
 			* 3 = Neitherside determiable (case C extended to F, or case E
			*/
			flag = 0;
			qblk = SGCB(syn,0,chrid,i);
            tblk = QT_SGCB(syn, qblk);
            tcon = QT_SGC(syn, qblk);
			if(missing){
				flag=3;
        		printf("%s\t%s\t.\t.\t%d\n",
               		seqname, tcon->name,flag);
					break;
			}


			qnode = cmap->map[qblk->linkid];
//			printf("[%d]\t%u::%u\t %u::%u\n",i,start,stop,qnode->feature->start,qnode->feature->stop);
			// Set return region assumes case D;
			tblk->start= qnode->match->start;
			tblk->stop = qnode->match->stop;
			// Start is BEFORE current query Block;
			if (start< qblk->start) {
				//Move down contiguous block, stoping at leftmost possible point
				while(qnode->prev != NULL && start > qnode->feature->stop){
					qnode = qnode->prev;
				}
				if(qnode->flag >1 && start < qnode->feature->stop){ //avoid duplicates in cases of overlap
					continue;
				}
				if(start < qnode->feature->start){ // Check we didn't advance into a case E,F Situation
					flag = 1;
					tblk ->start = qnode->match->start;
				} else if( start > qnode->feature->start)  {	//Start is contained within current block C,D
					tblk->start = qnode->match->start;
				} else {	//Case A,B situations
					tblk->start = qnode->match->stop;
				}
			
			}
			//Stop is AFTER current query Block
			if (stop > qblk->stop){ //Stop is AFTER
				while(qnode->next != NULL && stop > qnode->feature->start){
					qnode = qnode->next;
					//Advance i to avoid repeated ranges due to continuity.
					i++;
				}
				if(qnode->flag >1 && stop < qnode->feature->start){ //avoid duplicates in cases of overlap
					continue;
				}
				if(stop > qnode->feature->stop){ //Case E,F situations
					flag = flag==1? 3:2;
					tblk->stop = qnode->match->stop;					
				} else if (stop < qnode->feature->stop){
					tblk->stop = qnode->match->stop;
				} else {
					tblk->stop = qnode->match->start;
				}
			}
		
        	printf("%s\t%s\t%u\t%u\t%d\n",
               	seqname, tcon->name, tblk->start, tblk->stop, flag);
	
		}
        free(contigs->name);
        free(contigs->block);
        free(contigs);
	}
	
	free_contiguous_map(cmap);
}	

void free_contiguous_map(ContiguousMap * cmap){
	for(int i=0;i<cmap->size;i++){
		if(cmap->map[i]) free(cmap->map[i]);
	}
	if(cmap->map) free(cmap->map);
	if(cmap) free(cmap);
}


ContiguousMap * populate_contiguous_map(Synmap * syn){
	// Figure out how many unique blocks the current synmap has	
	size_t size = 0;
	for(uint32_t i=0; i< syn->genome[0]->size; i++){
    	size += syn->genome[0]->contig[i]->size;
	}
	
    // Initialize new ContiguousMap to that size to serve as hashmap
    // when looking up overlaping blocks
	ContiguousMap * cmap = init_contiguous_map(size);
	// Initialize ContiguousList structure to hold head and tail of current
	// Contiguous Set
	ContiguousList * clist = calloc(1,sizeof(ContiguousList));
	cmap->size = size;

	for(uint32_t i =0; i < syn->genome[0]->size; i++){
		
		//Grab contig to variable, to save writing
		Contig * ctig = (Contig *)malloc(sizeof(Contig));
		ctig = syn->genome[0]->contig[i];
		uint32_t forward_bound[2]= {0,0};	
		for(uint32_t j=0; j < ctig->size; j++){

			//set current block up to be a new node
			ContiguousNode *cnode = (ContiguousNode*)malloc(sizeof(ContiguousNode));
			cnode->next = NULL;
			cnode->prev = NULL;
			cnode->feature =ctig->block[j];
			cnode->match = syn->genome[1]->contig[cnode->feature->oseqid]->block[cnode->feature->oblkid];
			cmap->map[cnode->feature->linkid] = cnode;
			//clist first and last are null, as should only be the case during the first iteration of loop
			if(clist->first == NULL){
				contiguous_list_reset(clist,cnode);
				cmap->map[ctig->block[j]->linkid] = cnode;
//				printf("head \t%u::%u \t[%u:%u] \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid);
			//Otherwise we get to figure out if the current block can be added to or breaks the current ContiguousSet
			} else {
		
			//// Test for overlap on query side
				bool q_overlap = false;
				bool t_overlap = false;
				q_overlap = block_overlap(cnode->feature,clist->last->feature);
				t_overlap = (block_overlap(cnode->match,clist->last->match) && 
							 (cnode->feature->oseqid == clist->last->feature->oseqid));
				
				if(q_overlap){
					 cnode->flag=2;
//					 printf("query overlap \t%u::%u || %u::%u \t[%u:%u] \n",cnode->feature->start,cnode->feature->stop,
//						clist->last->feature->start, clist->last->feature->stop, cnode->feature->oseqid,cnode->feature->oblkid);
				}
				if(t_overlap){
					cnode->flag=3;
//					printf("target overlap \t%u::%u || %u::%u \t[%u:%u] \n",cnode->match->start,cnode->match->stop,
//						clist->last->match->start, clist->last->match->stop, cnode->feature->oseqid,cnode->feature->oblkid);
				}
				if(q_overlap || t_overlap){
					cmap->map[cnode->feature->linkid] = cnode;
					continue;
				}
			}
		// Setup integers to do comparison testing of targets
				
			int p_start = block_cmp_start(cnode->match, clist->last->match);
			int p_stop	= block_cmp_stop(cnode->match, clist->last->match);
			unsigned int p_t_ctig = clist->last->feature->oseqid;
					
			int n_start = 1;
			int n_stop	= 1;
			unsigned int n_t_ctig = cnode->feature->oseqid;

			if(j+1 < ctig->size){
				Block * nextBlock = (Block*)malloc(sizeof(Block));
				nextBlock = syn->genome[1]->contig[ctig->block[j+1]->oseqid]->block[ctig->block[j+1]->oblkid];
				n_start = block_cmp_start(nextBlock,cnode->match);
				n_stop = block_cmp_stop(nextBlock,cnode->match);
				n_t_ctig = ctig->block[j+1]->oseqid;
			}
			int p_diff_contig = (p_t_ctig == cnode->feature->oseqid)?false:true;	
			int n_diff_contig = (n_t_ctig == cnode->feature->oseqid)?false:true;;
			//Current feature target contig is not on same contig as previous target.
			if(p_diff_contig){
				contiguous_list_reset(clist,cnode);
				cmap->map[ctig->block[j]->linkid] = cnode;
			}
			// Test if the current block already exists in the cmap if it does, see if
			// there is flag indicating that it is immediately before or after an upcontig
			// twist
			if((j == forward_bound[0] || j == forward_bound[1]) && j!=0 ){
				cmap->map[ctig->block[j]->linkid] = cnode;

				if(j== forward_bound[1]){
					forward_bound[0]= 0;
					forward_bound[1]= 0;
				}
//				printf("uptwist breakpoint \t%u::%u \t[%u:%u] \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid);
			}
			// Test for Down Contig Twist
			if(p_start == -1 && !n_diff_contig){
				cnode->flag = -1;
				cmap->map[ctig->block[j]->linkid] = cnode;
				bool overlap = block_overlap(cmap->map[ctig->block[j]->linkid]->match,cmap->map[ctig->block[j]->linkid]->match);
				if(cmap->map[ctig->block[j-1]->linkid]->flag ==-1 && !overlap){
					//cnode->prev = cmap->map[ctig->block[j-1]->linkid];
					//cmap->map[ctig->block[j-1]->linkid]->next = cnode;;
					contiguous_list_push(clist,cnode);
					cmap->map[ctig->block[j-1]->linkid]->next = cnode;
					cmap->map[ctig->block[j]->linkid]->prev = cnode;
				}else{
					contiguous_list_reset(clist,cnode);
				}
//				printf("downtwist \t%u::%u \t[%u:%u] \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid);
			// Test for Next Contig Twist
			} else if ((n_start == -1 && n_stop == -1) && !n_diff_contig){
				int target_start=j;
				int target_current = target_start;
				int found =0;	
				// Find if next blocks are contiguous with query block
				Block * current = (Block*)malloc(sizeof(Block));
				Block * next = (Block*)malloc(sizeof(Block));
				while(j < ctig->size-1 && !found){
					current = syn->genome[1]->contig[ctig->block[j]->oseqid]->block[ctig->block[j]->oblkid];
					next = syn->genome[1]->contig[ctig->block[j+1]->oseqid]->block[ctig->block[j+1]->oblkid];
					//Contiguous set, rev	
					if  (block_cmp_start(next,current) == 1 && !block_overlap(next,current) 
						 && (ctig->block[j]->oseqid == ctig->block[j+1]->oseqid)
					    ){
					
						j++;
					} else {
						found = 1;
						target_current = target_start + 1;
						forward_bound[0]= target_current;
						
					}
				}
				
							
				current = clist->last->match;
				found =0;

				// Find where target is finally after the advanced mark 
				while(target_current < ctig->size -1 && found == 0){
					next = syn->genome[1]->contig[ctig->block[target_current+1]->oseqid]->block[ctig->block[target_current+1]->oblkid];
					if(block_cmp_start(next,current) !=1){
						target_current++;
					} else {
						found =1;
					 }
				}
				forward_bound[1] =target_current;	
				j= target_start;
				cmap->map[ctig->block[j]->linkid] = cnode;

//				printf("uptwist \t%u::%u \t[%u:%u] \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid);
			// Test for regular memeber of contiguous set
			}else if((p_start ==1 && p_stop == 1 && (n_stop ==1||n_start==1)) 
						|| (p_start==0 && p_stop ==0 && (j==0 || p_diff_contig))
					 	|| (p_start==1 && p_stop==1 && (n_diff_contig || j==ctig->size-1))
					){
			
				contiguous_list_push(clist,cnode);
				if(j>0){
					cmap->map[ctig->block[j-1]->linkid]->next = cnode;
					cmap->map[ctig->block[j]->linkid]->prev = cnode;
				}
//			printf("continous \t%u::%u \t[%u:%u] \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid);
			// Singletoen
			} else {
				contiguous_list_reset(clist,cnode);
				cmap->map[ctig->block[j]->linkid] = cnode;
//				printf("other \t%u::%u \t[%u:%u] \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid);
			}
		
		}
		clist->first = NULL;
		clist->last = NULL;
	}		
	//The actual list has done its job, free.
	free(clist);
	return cmap;

}

void contiguous_list_push(ContiguousList *clist, ContiguousNode *cnode){
	clist->last->next = cnode;
	clist->last = cnode;
}

void contiguous_list_reset(ContiguousList *clist, ContiguousNode *cnode){
	clist->first = cnode;
	clist->last = cnode;
}




