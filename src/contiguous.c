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
				if(i==0){
					region[0] = endblk->oblkid;
					region[1] = endblk->oblkid;
				} else {
					if(endblk->oblkid < region[0]){region[0]=endblk->oblkid;}
					if(endblk->oblkid > region[1]){region[1]=endblk->oblkid;}
				}
//        	printf("[%d]\t%u\t%u\t%d\n",
//              	qblk->linkid,qblk->start,qblk->stop, cmap->map[tblk->linkid]->flag);
//        	printf("[%d]\t%u\t%u\n",
//              	qblk->linkid,tblk->start,tblk->stop);
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
				if(qnode->flag >1 && start > qnode->feature->stop){ //avoid duplicates in cases of overlap
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
	cmap->size = size;
	
	for(uint32_t i =0; i < syn->genome[0]->size; i++){
		
		//Grab contig to variable, to save writing
		Contig * ctig = (Contig *)malloc(sizeof(Contig));
		ctig = syn->genome[0]->contig[i];
		uint32_t overlap_bound[2]= {0,0};	
		
		for(uint32_t j=0; j < ctig->size; j++){

			//set current block up to be a new node
			ContiguousNode *cnode = (ContiguousNode*)malloc(sizeof(ContiguousNode));
			cnode->flag =0;
			cnode->next = NULL;
			cnode->prev = NULL;
			cnode->feature =ctig->block[j];
			cnode->match = syn->genome[1]->contig[cnode->feature->oseqid]->block[cnode->feature->oblkid];
			//Start of contig is always a default node
				cmap->map[cnode->feature->linkid] = cnode;
			if(j==0){
				cmap->map[cnode->feature->linkid] = cnode;
				continue;
//				printf("head \t%u::%u \t[%u:%u] \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid);
			//Otherwise we get to figure out if the current block can be added to or breaks the current ContiguousSet
			} else {
		
			//// Test for overlap on query side
				bool q_overlap = false;
				bool t_overlap = false;
				for(int k = overlap_bound[0]; k<= overlap_bound[1]; k++){
					q_overlap = block_overlap(cnode->feature,ctig->block[k]) || q_overlap;
					t_overlap = (block_overlap(cnode->match,cmap->map[ctig->block[k]->linkid]->match) && 
							 (cnode->feature->oseqid == ctig->block[k]->oseqid)) || t_overlap;
					if(q_overlap || t_overlap) continue;
				}
					if(q_overlap){
						 cnode->flag=2;
//					 printf("[%u]\tquery overlap \n",j);
					}
					if(t_overlap){
						cnode->flag=3;
//					 	printf("[%u]\ttarget overlap \n",j);
					}
				if(q_overlap || t_overlap){
					if(q_overlap && t_overlap) cnode->flag = 4;
					overlap_bound[1]=j;
					cmap->map[cnode->feature->linkid] = cnode;
					continue;
				} else {
					overlap_bound[0] = j;
					overlap_bound[1] = j;
				}
			}
			// On different contig from previous 
			if(cnode->feature->oseqid != cmap->map[ctig->block[j-1]->linkid]->feature->oseqid){
//				printf("[%u]\tDifferent contig\n",j);
				cmap->map[cnode->feature->linkid] = cnode;
			} else if (cnode->feature->oblkid == ctig->block[j-1]->oblkid+1){ // Regular contiguous interval
				cmap->map[cnode->feature->linkid] = cnode;
				cmap->map[ctig->block[j-1]->linkid]->next = cnode;
				cmap->map[cnode->feature->linkid]->prev = cmap->map[ctig->block[j-1]->linkid];
//				printf("[%u]\tcontiguous\n",j);
			} else if( cnode->feature->oblkid < ctig->block[j-1]->oblkid){ //Twist to left of previous block
				cnode->flag = cmap->map[ctig->block[j-1]->linkid]->flag  < -1 ? -3:-1;
				cmap->map[cnode->feature->linkid] = cnode;
				if(cnode->feature->oblkid == ctig->block[j-1]->oblkid - 1){
					cmap->map[ctig->block[j-1]->linkid]->next = cnode;
					cmap->map[cnode->feature->linkid]->prev = cmap->map[ctig->block[j-1]->linkid];
				}
//					printf("[%u]\tTwisted Left\n",j);
			} else if( cnode->feature->oblkid > ctig->block[j-1]->oblkid){// Twist to right, possible transposition
				cnode->flag=-2;
				cmap->map[cnode->feature->linkid] = cnode;
//				printf("[%u]\tTwisted Right\n",j);
			} else { // Default case that should never be reached.	
				cmap->map[cnode->feature->linkid] = cnode;
//				printf("[%u]\tOTHER\n",j);
			}	
		}

	}		
	//The actual list has done its job, free.
	return cmap;

}




