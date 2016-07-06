#ifndef __LEV_H__
#define __LEV_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#define TOLC(s)	for(char *p =(s);*p;++p) *p=*p>0x40&&*p<0x5b?*p|0x60:*p
#define LEVMIN(x,y,z) ((x)<(y)?((x)<(z)?(z):(z)):((y)<(z)?(y):(z)))

/** 
 * @brief Dictionary structure for contig names
 *
 */
typedef struct DictNode {
  char *word;   /**< The  word. */
  uint32_t arr_loc; /**< Word's mapped index in the synteny db */
  bool original; /**< If this word was taken from synteny db */
  struct DictNode *next; /**< Next dictionary entry */
} DictNode;

/**
 * @brief Calculates levenshtein distance between two strings
 *
 * This function calculates the levenshtein disance between two strings
 * using an unweighed, simple dynamic programming approach. It is used
 * to see if a given string matches a member in the dictionary, or to provide
 * the closest match.
 *
 * @param char* s1 First string
 * @param char* s2 Second string
 *
 * @return int32_t distance distance between two strings
 */
int32_t lev_dist(char *s1, char *s2);

/**
 * @brief Calculates levenshtein distance between two strings
 *
 * This function calculates the levenshtein disance between two strings
 * using an unweighed, simple dynamic programming approach. It is used
 * to see if a given string matches a member in the dictionary, or to provide
 * the closest match.
 *
 * @param FILE* synfile	Synteny database file as prepared  by synder
 * @param FILE* intfile Gff file to be converted
 * @param int swap Switch to indicate if matching names on target or query side
 *
 */
void convert_seqname(FILE * synfile, FILE * intfile, int swap);

#endif
