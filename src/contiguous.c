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

void contiguous_query(Synmap * syn, FILE * intfile, bool pblock){
	// Ensure blocks are sorted (BAD IDEA, REDACTED FOR NOW)
	//	sort_all_contigs(syn);
	// count total number of unique block-block pairs for hashmap
	ContiguousMap *cmap= populate_contiguous_map(syn);
    
	uint32_t interval =0;
	uint32_t missloc =0;
	char seqname[128];
    int chrid, start, stop,flag;
    Contig * contigs;
    Contig * tcon;
    Block * qblk;
	uint32_t blkid;
	ContiguousNode * qnode;
	ContiguousNode * original;
	Block twoblk;
    bool missing;
    size_t length = 1024;
	char *line = (char*) malloc(length*sizeof(char));

while(fgets(line,length,intfile) && !feof(intfile)){
	if(!sscanf(line, "%d %*s %*s %d %d %*s %*c %*s %s\n", &chrid, &start, &stop, seqname)){
		printf("invalid input\n");
		continue;
	}
        contigs = get_region(SGC(syn, 0, chrid), start, stop);
		uint32_t region[2]={0,0};
        missing = false;

	    if(contigs->size == 2){
 	      twoblk.start = start;
	        twoblk.stop = stop;
    	    if(!((CB(contigs, 0) && block_overlap(CB(contigs, 0), &twoblk)) ||
        	     (CB(contigs, 1) && block_overlap(CB(contigs, 1), &twoblk))
            	)){
					if(contigs->block[0]){
						missloc= cmap->map[contigs->block[0]->linkid]->qblkid;
					} else {
						missloc= cmap->map[contigs->block[1]->linkid]->qblkid;
						missloc = missloc >0 ? missloc : 0;

            		}
					missing = true;
        	}
     	}

        for(int i = 0; i < contigs->size; i++){
            qblk = contigs->block[i];
            if(qblk){
				blkid = cmap->map[qblk->linkid]->qblkid;
				if(i==0){
					region[0] = blkid;
					region[1] = blkid;
				} else {
					if(blkid < region[0]){region[0]=blkid;}
					if(blkid > region[1]){region[1]=blkid;}
				}
//        	printf("[%d]\t%u\t%u\t\n",
//              	qblk->linkid,qblk->start,qblk->stop);
//        	printf("[%d]\t{%u}\n",
//              	qblk->linkid,blkid);
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
			Contig* qcon = SGC(syn,0,chrid);
            Block* tblk = init_block( QT_SGCB(syn, qblk)->start, QT_SGCB(syn,qblk)->stop,0,0,0);
            tcon = QT_SGC(syn, qblk);
			Block* q_blk;
			Block* t_blk;
			Contig* t_con;
			
			if(missing){
				if(pblock){		
					q_blk = SGCB(syn,0,chrid,missloc);
					t_blk = QT_SGCB(syn,q_blk);
					t_con = QT_SGC(syn,q_blk);
      				printf("Q\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\n",
            		   	seqname, qcon->name, q_blk->start, q_blk->stop,
						t_con->name,t_blk->start,t_blk->stop, interval);
					q_blk = SGCB(syn,0,chrid,missloc+1);
					t_blk = QT_SGCB(syn,q_blk);
					t_con = QT_SGC(syn,q_blk);
      				printf("Q\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\n",
            		   	seqname, qcon->name, q_blk->start, q_blk->stop,
						t_con->name,t_blk->start,t_blk->stop, interval);

				}
				q_blk = SGCB(syn,0,chrid,missloc);
				t_blk = QT_SGCB(syn,q_blk);
				t_con = QT_SGC(syn,q_blk);
				int64_t offset;
				if(start < q_blk->start){ // query region is before block
					if(cmap->map[q_blk->linkid]->flag >-2){ 
					// return from start of block, to offest to start of query
					// on target side
						flag = 4;
					    tblk->stop = t_blk->start;
						offset = t_blk->start - (q_blk->start-start);
					    tblk->start = offset < tblk->stop && offset > 0 ? (uint32_t)offset : 0;
					} else {
						flag = 10;
						tblk->start = t_blk->stop; 
						tblk->stop = t_blk->stop + (q_blk->start - start);
					}
				} else { // query region after block
					
					if(cmap->map[q_blk->linkid]->flag >-2){
					// return from end of block, to offest to end of query
					// on target side
						flag = 11;
						tblk->start = t_blk->stop; 
						tblk->stop = t_blk->stop + (stop - q_blk->stop);
					} else {
						flag = 4;
					    tblk->stop = t_blk->start;
						offset = t_blk->start - (stop-q_blk->stop);
					    tblk->start = offset < tblk->stop && offset > 0 ? (uint32_t)offset : 0;
					}
				}
				
				printf(">\t%u\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%d\n",
            	   			interval,seqname,qcon->name,start,stop,
							tcon->name, tblk->start,tblk->stop,flag);
				
				blkid = cmap->map[q_blk->linkid]-> qblkid;
				if (start < q_blk->start && blkid >0){
					q_blk = SGCB(syn,0,chrid,blkid-1);
					t_blk = QT_SGCB(syn,q_blk);
            		tcon = QT_SGC(syn, q_blk);

					if(cmap->map[q_blk->linkid]->flag >-2){
						flag = 12;
					    tblk->start = t_blk->stop;
						offset = t_blk->stop + (stop-q_blk->stop);
					    tblk->stop = (uint32_t)offset>=0 ?offset:0;
					} else {
						flag = 4;
					    tblk->stop = t_blk->start;
						offset = t_blk->start - (stop-q_blk->stop);
					    tblk->start = offset < tblk->stop && offset > 0 ? (uint32_t)offset : 0;
					}	
			    	printf(">\t%u\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%d\n",
        	   			interval,seqname,qcon->name,start,stop,
						tcon->name, tblk->start,tblk->stop,flag);
				
				} else if(blkid +1 < qcon->size){ // query region after block
					q_blk = SGCB(syn,0,chrid,blkid+1);
					t_blk = QT_SGCB(syn,q_blk);
            		tcon = QT_SGC(syn, q_blk);
					
					if(cmap->map[q_blk->linkid]->flag >-2){
						flag = 4;
					    tblk->stop = t_blk->start;
						offset = t_blk->start - (q_blk->start-start);
					    tblk->start = offset < tblk->stop && offset > 0 ? (uint32_t)offset : 0;
					} else {
						flag = 13;
					    tblk->start = t_blk->stop;
						offset = t_blk->stop + (q_blk->start - start);
					    tblk->stop = (uint32_t)offset>=0 ?offset:0;
					}	
					
			    	printf(">\t%u\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%d\n",
        	   			interval,seqname,qcon->name,start,stop,
						tcon->name, tblk->start,tblk->stop,flag);
				}
			
    // query region is before block

				interval++;
				free_block(tblk);
				break;
			}
			
			qnode = cmap->map[qblk->linkid];
		    original = cmap->map[qblk->linkid];

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
					if(qnode->flag > -2  || (original->next != NULL && original->next->flag == 0)){ 
						// The next check here and in the next block is to detect edge cases where block we are checking is
						// registering as a left translation, but is part of a continuous run to the right
						// most often occcurs in the middle of messy overlap blocks

						flag = 1;
						tblk ->start = qnode->match->start;
					} else {
						flag = 2;
						tblk ->stop = qnode->match->stop;
					}
				} else if( start > qnode->feature->start)  {	//Start is contained within current block C,D
					if(qnode->flag > -2){
						tblk->start = qnode->match->start;
					} else {
						tblk ->stop = qnode->match->stop;
					}
				} else {	//Case A,B situations
					if(qnode->flag > -2){
						tblk->start = qnode->match->stop;
					} else {
						tblk ->stop = qnode->match->start;
					}
				}
			
			}
			if(pblock){	
				//qnode = cmap->map[qblk->linkid];
				q_blk = SGCB(syn,0,chrid,qnode->qblkid >0?qnode->qblkid-1:0);
				t_blk = QT_SGCB(syn,q_blk);
				t_con = QT_SGC(syn,q_blk);
      			printf("Q\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\t\n",
            	   	seqname, qcon->name, q_blk->start, q_blk->stop,
					t_con->name,t_blk->start,t_blk->stop, interval);
				t_blk = SGCB(syn,1,qnode->feature->oseqid,
					(qnode->feature->oblkid >0 ? qnode->feature->oblkid-1:0));
				q_blk = SGCB(syn,0,t_blk->oseqid, t_blk->oblkid);
				t_con = SGC(syn,1,q_blk->oseqid);
      			printf("T\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\n",
            	   	seqname, qcon->name, q_blk->start, q_blk->stop,
					t_con->name,t_blk->start,t_blk->stop, interval);
      	
      			printf("I\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\n",
            	   	seqname, qcon->name, qblk->start, qblk->stop,
					tcon->name,qnode->match->start,qnode->match->stop, interval);
			}			

			
			//Stop is AFTER current query Block
			if (stop > qblk->stop){ //Stop is AFTER
				while(qnode->next != NULL && stop){
					ContiguousNode* temp = qnode->next;
					//Advance i to avoid repeated ranges due to continuity.
					if(stop > temp->feature->start){	
						i++;
						qnode = temp;
						if(pblock){
							t_con = SGC(syn,1,qnode->feature->oseqid);
      						printf("I\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\n",
               					seqname, qcon->name, qnode->feature->start, qnode->feature->stop,
								t_con->name,qnode->match->start,qnode->match->stop, interval);
						}
					} else {
						break;
					}
				}
				if(qnode->flag >1 && stop < qnode->feature->start){ //avoid duplicates in cases of overlap
					continue;
				}
				if(stop > qnode->feature->stop){
					if(qnode->next == NULL){ //Case E,F situations
						if(qnode->flag > -2){
							flag = flag==1? 3:2;
							tblk->stop = qnode->match->stop;					
						} else {
							flag = flag==2? 3:1;
							tblk ->start = qnode->match->start;
						}
					} else {	//Case A,B situations
						q_blk = SGCB(syn,0,chrid, i+1 < qcon->size ? i+1:i);
						t_blk = QT_SGCB(syn,q_blk);
						if(cmap->map[qblk->linkid]->flag > -2 || qnode->next->flag == 0){ //similar check to above 
							tblk->stop = t_blk->start;
							
						} else {
							tblk->start = t_blk->stop;
						}
					}
				} else { //Case C,D 
					tblk->stop = qnode->match->stop;
					if(qnode->flag > -2){
						tblk ->stop = qnode->match->stop;
					} else {
						tblk->start = qnode->match->stop;
					}
				}	

			}
			if(pblock){	
				if(stop < qnode->feature->start){	
					q_blk = qnode->feature;
				} else {
					q_blk = SGCB(syn,0,chrid, i+1 < qcon->size ? i+1:i);
				}
		
				t_blk = QT_SGCB(syn,q_blk);
				t_con = QT_SGC(syn,q_blk);
	      		printf("Q\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\n",
	               	seqname, qcon->name, q_blk->start, q_blk->stop,
					t_con->name,t_blk->start,t_blk->stop, interval);
				t_blk = SGCB(syn,1,qnode->feature->oseqid,qnode->feature->oblkid+1< tcon->size? qnode->feature->oblkid+1:qnode->feature->oblkid);
				q_blk = SGCB(syn,0,t_blk->oseqid, t_blk->oblkid);
				t_con = SGC(syn,1,q_blk->oseqid);
	      		printf("T\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%u\n",
	               	seqname, qcon->name, q_blk->start, q_blk->stop,
					t_con->name,t_blk->start,t_blk->stop, interval);
      		}
			printf(">\t%u\t%s\t%s\t%u\t%u\t%s\t%u\t%u\t%d\n",
               		interval,seqname,qcon->name,start,stop,
					tcon->name, tblk->start,tblk->stop,flag);

			interval++;
	
			free_block(tblk);
		}
        free(contigs->name);
        free(contigs->block);
        free(contigs);
	}
	free(line);
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
			cnode->qblkid = j;
			//Start of contig is always a default node
				cmap->map[cnode->feature->linkid] = cnode;
			if(j==0){
				cmap->map[cnode->feature->linkid] = cnode;
//				printf("%d\t head \t%u::%u \t[%u:%u] {%u,%u} \n",i,cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid,j,cnode->match->oblkid);
//				printf("\t[%u:%u]  \n",cnode->match->start,cnode->match->stop);
				continue;
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
//						printf("Query Overlap \t%u::%u \t[%u:%u] {%u,%u} \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid,j,cnode->match->oblkid);
//						printf("\t[%u:%u]  \n",cnode->match->start,cnode->match->stop);
					}
					if(t_overlap){
						cnode->flag=3;
//						printf("Target Overlap \t%u::%u \t[%u:%u] {%u,%u} \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid,j,cnode->match->oblkid);
//						printf("\t[%u:%u]  \n",cnode->match->start,cnode->match->stop);
//						printf("[Target] \t%u::%u \n",cnode->match->start,cnode->match->stop);
					}
				if(q_overlap || t_overlap){
					if(q_overlap && t_overlap) cnode->flag = 4;
					//if(cnode->feature->oblkid < ctig->block[j-1]->oblkid) cnode->flag = -2;
					overlap_bound[1]=j;
					cmap->map[cnode->feature->linkid] = cnode;
					continue;
				} 
				 else {
					overlap_bound[0] = j;
					overlap_bound[1] = j;
				}
			}
			// On different contig from previous 
			if(cnode->feature->oseqid != cmap->map[ctig->block[j-1]->linkid]->feature->oseqid){
//				printf("Different Contig \t%u::%u \t[%u:%u] {%u,%u} \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid,j,cnode->match->oblkid);
//				printf("\t[%u:%u]  \n",cnode->match->start,cnode->match->stop);
				cmap->map[cnode->feature->linkid] = cnode;
			} else if (cnode->feature->oblkid == ctig->block[j-1]->oblkid+1){ // Regular contiguous interval
				cnode->flag = 0;
				cmap->map[cnode->feature->linkid] = cnode;
				cmap->map[ctig->block[j-1]->linkid]->next = cnode;
				cmap->map[cnode->feature->linkid]->prev = cmap->map[ctig->block[j-1]->linkid];
//				printf("Contiguous \t%u::%u \t[%u:%u] {%u,%u} \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid,j,cnode->match->oblkid);
//				printf("\t[%u:%u]  \n",cnode->match->start,cnode->match->stop);
			} else if( cnode->feature->oblkid < ctig->block[j-1]->oblkid){ //Twist to left of previous block
				cnode->flag = (cmap->map[ctig->block[j-1]->linkid]->flag  == -2 || cmap->map[ctig->block[j-1]->linkid]->flag == 0) ? -2:0;
				cmap->map[ctig->block[j-1]->linkid]->next = NULL;
				cmap->map[cnode->feature->linkid] = cnode;
				if(cnode->feature->oblkid == ctig->block[j-1]->oblkid - 1){
					cnode->flag = -3;
					cmap->map[ctig->block[j-1]->linkid]->next = cnode;
					cmap->map[ctig->block[j-1]->linkid]->flag = -3;
					cmap->map[cnode->feature->linkid]->prev = cmap->map[ctig->block[j-1]->linkid];
				}
//				printf("Twist Left \t%u::%u \t[%u:%u] {%u,%u} \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid,j,cnode->match->oblkid);
//				printf("\t[%u:%u]  \n",cnode->match->start,cnode->match->stop);
			} else if( cnode->feature->oblkid > ctig->block[j-1]->oblkid){// Twist to right, possible transposition
				if(j+1< ctig->size){
				    cnode->flag= ctig->block[j]->oblkid == ctig->block[j+1]->oblkid+1 ? -3:-1;
				} else {
					cnode->flag = -1;
				}
				cmap->map[cnode->feature->linkid] = cnode;
//				printf("Twist Right \t%u::%u \t[%u:%u] {%u,%u} \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid,j,cnode->match->oblkid);
//				printf("\t[%u:%u]  \n",cnode->match->start,cnode->match->stop);
			} else { // Default case that should never be reached.	
				cmap->map[cnode->feature->linkid] = cnode;
//				printf("OTHER \t%u::%u \t[%u:%u] {%u,%u} \n",cnode->feature->start,cnode->feature->stop,cnode->feature->oseqid,cnode->feature->oblkid,j,cnode->match->oblkid);
//				printf("\t[%u:%u]  \n",cnode->match->start,cnode->match->stop);
			}	
		}

	}		
	//The actual list has done its job, free.
	return cmap;

}



