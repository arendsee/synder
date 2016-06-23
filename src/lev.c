#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "lev.h"

void convert_seqname(FILE* synfile, FILE* intfile){

    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char seqid[128];
    uint32_t ncontigs,id;
	uint32_t line_no= 0;
    char **name_arr;
   
    DictNode *root = (DictNode*)malloc(sizeof(DictNode));
    DictNode *current;
    current = root;
    
	while ((read = getline(&line, &len, synfile)) != EOF) {
        line_no++;
        if(line[0] == '>'){
            if(sscanf(line, "> %s %u %*c", seqid, &ncontigs)){
                name_arr = (char**)malloc(ncontigs*sizeof(char*));
			}
        }
        else if(line[0] == '@'){
            break;
        }
        else if(line[0] == '$') {
            if(sscanf(line, "$ %u %*u %s %*c\n", &id, seqid)){
                TOLC(seqid);
                name_arr[id] = (char*)malloc(strlen(seqid));
		        strcpy(name_arr[id],seqid);
			    current-> word = (char*)malloc(strlen(seqid));
			    strcpy(current->word,seqid);
			    current-> arr_loc = id;
				current->original= true;
                current->next = (DictNode*)malloc(sizeof(DictNode));
                current= current->next;
			}
				
        }
        else {
            fprintf(stderr, "Incorrect file format, line %d\n", line_no);
            exit(EXIT_FAILURE);
        }
    }
	free(current->next);
    current =root;

    char seqname[128];
    int start,end;
    char strand;
	char attribute[512];
    size_t length = 1024;
    char *nline = (char*) malloc(length*sizeof(char));

    while(fgets(nline,length,intfile) && !feof(intfile)){
        if(!sscanf(nline, "%s %*s %*s %d %d %*s %c %*s %s\n", seqname, &start, &end, &strand, attribute)){
            continue;
        }

        uint32_t dist;
        uint32_t min_idx=0;
        uint32_t idx = 0;
        uint32_t min=1000;
		TOLC(seqname);
    	current =root;
		while(current->next != NULL){
			idx=current->arr_loc;
			dist = lev_dist(seqname,current->word);
			if(dist == 0){
                min = dist;
            	min_idx = idx;
				break;
			} else if (dist < min && current->original){
            	min_idx = idx;
                min = dist;
            }

        	current= current->next;
		}
        if(min != 0){
			fprintf(stderr,"\n%s does not match any names, closest match: %s\n", seqname,name_arr[min_idx]);
			fprintf(stderr,"Press Enter to use the closest match, or enter a number from below: \n");
            for(int i = 0; i<= idx ;i++){
		    	fprintf(stderr,"%d: %s\t\t",i,name_arr[i]);
				if(i%5==4){
		    	    fprintf(stderr,"\n");
				}
			}
			
			uint32_t input = min_idx;
			char *ptr;
	     	while(true){
			    fprintf(stderr,"\n[enter] or [0-%d]: ",idx);
				fgets(line,30,stdin);
           	    if (line[0] == '\n'){
					min_idx = input;
					break;
				}
				input = strtoul(line,&ptr,10);
				
				if(input < idx){
					min_idx = input;
					break;
				}
				fprintf(stderr,"%d %u is an invalid choice, please Try again.\n",atoi(line),input);
			}

		}

		DictNode* newRoot = (DictNode*)malloc(sizeof(DictNode));
		newRoot->word  = (char*)malloc(strlen(seqname));
		strcpy(newRoot->word,seqname);
		newRoot-> arr_loc = min_idx;
		newRoot->original = false;
		newRoot->next = root;
		root = newRoot;
	
		printf("%u\t.\t.\t%d\t%d\t.\t%c\t.\t%s\n",min_idx, start, end, strand, attribute);
		
	}
}

int32_t lev_dist(char *s1, char *s2){

    uint32_t s1len, s2len;
    uint32_t curr,prev;

    s1len= strlen(s1);
    s2len= strlen(s2);

    uint32_t loc[s1len+1];

    //initialize
    uint32_t i,j;
    for(i = 1; i <= s1len; i++){
        loc[i]=i;
    }

    //horray for Dynamic Programming!
    for(i = 1; i <= s2len; i++){
        loc[0]= i;
        for(j =1, curr = i-1 ; j<=s1len; j++){
            prev = loc[j];
            curr += (s1[j-1]== s2[i-1] ? 0:1);
            loc[j] = LEVMIN(loc[j]+1, loc[j-1]+1,curr);
            curr = prev;
        }
    }

    return(loc[s1len]);
}

