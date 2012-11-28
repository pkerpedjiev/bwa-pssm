// PSSMSearcher: A library for searching with Position Specific Scoring Matrices
//
// Copyright (C) 2006, 2007 Jes Frellsen, Ida Moltke and Martin Thiim
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "shared.h"
#include "pssm.h"
#include "bwtaln.h"

#define min(x, y) (x < y ? x : y)


// Function prototypes
inline int map(const unsigned char *letters, int length, int alphabet_size);
inline int powI(int x, int y);



// ******************************************************************
// Functions for creating and mutating a pssm
// ******************************************************************

// Function that creates and allocates a matrix (except for the scores)
PSSM init_matrix(int order, int length, int alphabet_size){

  // Check if longer than max-length of PSSMs
  if(length >= MAXPSSMSIZE) {
    fprintf(stderr,"Matrix is to long.");
    return NULL;
  }

  // All attributes are initialised
  PSSM pssm = malloc(sizeof(struct PSSM));
  if(!pssm){fprintf(stderr,"Couldn't allocate memory for the PSSM.");return NULL;};
  memset(pssm, 0, sizeof(struct PSSM));

  pssm->order = order;
  pssm->length = length;
  pssm->alphabet_size = alphabet_size;
  
  pssm->thresholds = malloc((length)*sizeof(int));
  if(!pssm->thresholds) {free(pssm);fprintf(stderr,"Couldn't allocate memory for the thresholds.");return NULL;}

  pssm->bi = malloc((length)*sizeof(int));
  if(!pssm->bi) {free(pssm);fprintf(stderr,"Couldn't allocate memory for the best inexact scores.");return NULL;}
 
  pssm->be = malloc((length)*sizeof(int));
  if(!pssm->be) {free(pssm);fprintf(stderr,"Couldn't allocate memory for the best exact scores.");return NULL;}
 

  pssm->saved_scores = malloc((length+1)*sizeof(int));
  if(!pssm->saved_scores) {free(pssm);fprintf(stderr,"Couldn't allocate memory for the saved_scores.");return NULL;}
 
  pssm->offsets = malloc((length+1)*sizeof(int));
  if(!pssm->offsets) {free(pssm);fprintf(stderr,"Couldn't allocate memory for the offsets.");return NULL;}

  calc_and_set_offsets(pssm,order,length,alphabet_size);

  pssm->scores = malloc(pssm->offsets[length]*sizeof(int));
  if(!pssm->scores) {free(pssm);fprintf(stderr,"Couldn't allocate memory for the PSSM.");return NULL;}
  memset(pssm->scores, 0, pssm->offsets[length]*sizeof(int));
  memset(pssm->saved_scores,0,(length+1) * sizeof(int));


  // Allocate space for the scores
  //fprintf(stderr, "length:%d pssm->offsets[length]:%d\n", length, pssm->offsets[length]);
  //memset(pssm->saved_scores, 0, MAXPSSMSIZE+1 * sizeof(int));

  return pssm;
}


// Creates a matrix given an array of scores.
//
// The function also checks if the length is to high and if numbers of
// scores in '*scores' array matches order, length and alphabet_size
// (returns NULL is check fails)
PSSM init_matrix_score(int order, int length, int alphabet_size, int *scores, int nScores, int threshold){

  // Check if longer than max-length of PSSMs
  if(length >= MAXPSSMSIZE) {
    fprintf(stderr,"Matrix is to long.");
    return NULL;
  }

  // All attributes are initialised
  PSSM pssm = malloc(sizeof(struct PSSM));
  memset(pssm, 0, sizeof(struct PSSM));

  if(!pssm){fprintf(stderr,"Couldn't allocate memory for the PSSM.");return NULL;};
  pssm->order = order;
  pssm->length = length;
  pssm->alphabet_size = alphabet_size;
  pssm->scores = scores;
 
  pssm->thresholds = malloc((length)*sizeof(int));
  if(!pssm->thresholds) {free(pssm);fprintf(stderr,"Couldn't allocate memory for the thresholds.");return NULL;}
 
  pssm->bi = malloc((length)*sizeof(int));
  if(!pssm->bi) {free(pssm);fprintf(stderr,"Couldn't allocate memory for the best inexact scores.");return NULL;}
 
  pssm->be = malloc((length)*sizeof(int));
  if(!pssm->be) {free(pssm);fprintf(stderr,"Couldn't allocate memory for the best exact scores.");return NULL;}


  pssm->saved_scores = malloc((length+1)*sizeof(int));
  if(!pssm->saved_scores) {free(pssm);fprintf(stderr,"Couldn't allocate memory for the saved_scores.");return NULL;}
 
  pssm->offsets = malloc((length+1)*sizeof(int));
  if(!pssm->offsets) {free(pssm);fprintf(stderr,"Couldn't allocate memory for the offsets.");return NULL;}


  memset(pssm->saved_scores,0,(length+1) * sizeof(int));

  calc_and_set_offsets(pssm,order,length,alphabet_size);

  // Check if length of scores is correct
  if(nScores != pssm->offsets[length]) {
    free(pssm);
    char errormsg[160];
    sprintf(errormsg, "Mismatch between list size (%i) and size calculated from order, length and alphabet size (%i).",
	    nScores, pssm->offsets[length]);
    fprintf(stderr,errormsg);
    return NULL;
  }

  calc_and_set_thresholds(pssm, threshold);

  return pssm;
}


// Function that frees the memory of a matrix
void release_matrix(PSSM pssm){
  if(pssm){
    if(pssm->scores)
      free(pssm->scores);
    if (pssm->offsets)
        free(pssm->offsets);
    if (pssm->thresholds)
        free(pssm->thresholds);
    if (pssm->bi)
        free(pssm->bi);
    if (pssm->be)
        free(pssm->be);
    if (pssm->saved_scores)
        free(pssm->saved_scores);
    free(pssm);
  }
}


// Function that fills in a given pssms offset array
void calc_and_set_offsets(PSSM pssm, int order, int length, int alphabet_size){
  int i;
  int offset, realorder;

  // The offset for position 0 is set to 0
  offset = 0;
  pssm->offsets[offset] = 0;

  // The offset for position 'i' is set to the offset for posi-
  // tion 'i-1' plus the number of scores for position 'i-1'
  for(i = 1; i <= length; i++){
    realorder = min(i,order+1);
    offset += powI(alphabet_size,realorder);
    pssm->offsets[i] = offset;
  }
}


// Find max score for position pos
int maximum_score(PSSM pssm, int pos) {
	int key, offset, numScores;
	int curMax, curScore;

	// For each position 'pos' (starting with the last)
    curMax = -INT_MAX;
    numScores = pssm->offsets[pos+1] - pssm->offsets[pos];
    offset = pssm->offsets[pos];

    // The highest score is found
    for(key = 0; key < numScores; key++) {
    	curScore = pssm->scores[offset+key];
    	if(curScore > curMax) curMax = curScore;
    }
    return curMax;
}


// Function that fills in a given pssms threshold array
void calc_and_set_thresholds(PSSM pssm, int threshold){
	int pos;

	// For each position 'pos' (starting with the last)
    // The threshold for 'pos' is set to 'threshold' minus
    // the sum of the max scores on all positions after
	pos = get_length(pssm)-1;
	pssm->thresholds[pos] = threshold;
	for( --pos ; pos >= 0; pos--){
		pssm->thresholds[pos] = pssm->thresholds[pos+1]-maximum_score(pssm, pos+1);
	}

/*
    for (pos=0; pos < 21; pos++) {
       fprintf(stderr, "threshold[%d]: %f\n", pos, pssm->thresholds[pos]); 
    }
    */
}

void calculate_reverse_best_inexact(PSSM mat, int *min_drops)
{
    int min_drop = min_drops[0];
    int best_score = 0.0;
    int i;

    for (i = 0; i < mat->length; i++) {
        if (min_drops[i] < min_drop)
            min_drop = min_drops[i];
        best_score += maximum_score(mat, i);
        mat->bi[i] = best_score - min_drop;
    }
}

void calculate_reverse_best_exact(PSSM mat)
{
    int best_score = 0.0;
    int i;

    for (i = 0; i < mat->length; i++) {
        best_score += maximum_score(mat, i);
        mat->be[i] = best_score;
    }
}


// Function that fills in a given pssms threshold array
void calc_and_set_reverse_thresholds(PSSM pssm, int start, int end, int threshold){
	int i;

	// For each position 'pos' (starting with the last)
    // The threshold for 'pos' is set to 'threshold' minus
    // the sum of the max scores on all positions after
	//pos = get_length(pssm)-1;
	//pssm->thresholds[pos] = threshold;
    pssm->thresholds[start-1] = threshold;
	for( i = start ; i < end ; i++){
		//pssm->thresholds[pos] = pssm->thresholds[pos+1]-maximum_score(pssm, pos+1);
        pssm->thresholds[i] = pssm->thresholds[i-1] - maximum_score(pssm, i-1);
        //fprintf(stderr, "thresholds[%d] = %f\n", i, pssm->thresholds[i]);
	}

/*
    for (pos=0; pos < 21; pos++) {
       fprintf(stderr, "threshold[%d]: %f\n", pos, pssm->thresholds[pos]); 
    }
    */
}

void printThresholds(PSSM pssm)
{
    int i;
    for (i = 0; i < pssm->length; i++) {
        fprintf(stderr, "thresholds[%d] = %d\n", i, pssm->thresholds[i]);
    }
}

void print_horizontalPSSM(PSSM pssm, int n)
{
    int i;
    ubyte_t j;

    if (n == 0) 
        n = pssm->length;

    for (j = 0; j < 4; j++) {
        for (i = 0; i < n; i++)
            fprintf(stderr, "%6d ", (int)get_score_fast(pssm, &j, i));
        fprintf(stderr, "\n");
    }
}

void printPSSM(PSSM pssm)
{
    int i;
    ubyte_t j;
    for (i = 0; i < pssm->length; i++) {
        fprintf(stderr, "i: %d ", i);
        for (j = 0; j < 4; j++) {
            fprintf(stderr, " %df", get_score_fast(pssm, &j, i));
        }
        fprintf(stderr, "\n");
    }
}

// Function that sets a specific score in a given pssm
void set_score(PSSM pssm, const unsigned char *letters, int pos, int score){
  int index = 0;
  int realorder = min(pos, pssm->order);

  index =  pssm->offsets[pos];
  index += map(letters, realorder+1, pssm->alphabet_size);

  pssm->scores[index] = score;
}


// ******************************************************************
// Accessors for pssms
// ******************************************************************

// Function that returns the score of a given letter in a suffix
int get_score(PSSM pssm, const unsigned char *base_letter, int pos){
  int index = 0;
  int realorder = min(pos, pssm->order);
  unsigned char *letters = (unsigned char *)(base_letter-realorder);

  // The index of the letters are calculated
  index =  pssm->offsets[pos];
  index += map(letters, realorder+1, pssm->alphabet_size);

  //printf("pssm->score: %f", pssm->scores[index]);

  return pssm->scores[index];
}


// Function that returns the threshold of a given pos in a given pssm
int get_threshold(PSSM pssm, int pos){
  return pssm->thresholds[pos];
}


// Function that returns the global threshold for a given pssm
int get_global_threshold(PSSM pssm){
  return pssm->thresholds[pssm->length - 1];
}


// Function that returns the length of a given pssm
int get_length(PSSM pssm){
  return pssm->length;
}

// ******************************************************************
// Additional functions
// ******************************************************************

// Function for calculation small powers of intergers
int powI(int x, int y) {
  int i;
  int res = 1;

  for(i = 1; i<=y; i++)
    res = res*x;

  return res;
}


// Fuctions for mapping from char array to integer
inline int map(const unsigned char *letters, int length, int alphabet_size) {
  int k;
  int value = 0;

  for(k = length-1; k >= 0; k--){
    value += (*(letters+k))*powI(alphabet_size,(length-1)-k);
  }

  return value;
}

/* Read JASPAR-formatted PSSM from a file */
PSSM ReadPSSMFromFile(const char *filename)
{
    FILE *f = fopen(filename, "r");
    char *pch;
    char line[MAX_LINE_LENGTH];
    char alphabet[MAX_LINE_LENGTH];
    char *dAlphabet;
    int alphabet_size=0, pssmLength=0, counter=0;
    int *scores = (int *) malloc(MAX_LINE_LENGTH * MAXPSSMSIZE * sizeof(int));
    int *pssmScores;
    int i, j, nScores=0;
    int maxScore=0, maxRowScore=0;
    PSSM pssm;

    if (!f) {
        fprintf(stderr, "Read PSSM: Unable to open file: %s\n", filename);
        exit(1);
    }

    memset(alphabet, 0, MAX_LINE_LENGTH);

    while(fgets(line, MAX_LINE_LENGTH, f)) {
        if (alphabet_size >= MAX_LINE_LENGTH) {
            fprintf(stderr, "PSSM alphabet is too large. Maximum size is: %d\n", MAX_LINE_LENGTH);
            exit(1);
        }

       pch = strtok(line, " \t");
       
        alphabet[alphabet_size] = *pch;
        pch = strtok(NULL, " []\t");
        counter = 0;
        while (pch) {
            if (counter > MAXPSSMSIZE) {
                fprintf(stderr, "The PSSM is too long: Maximum length is: %d\n", MAXPSSMSIZE);
                exit(1);
            }
            scores[pssmLength * alphabet_size + counter] = atof(pch);

	    pch = strtok(NULL, " []\t");
            counter++;
        }
        if (pssmLength == 0)
            pssmLength = counter;

        alphabet_size++;
    }

    pssmScores = (int *)malloc(alphabet_size * pssmLength * sizeof(int));
    dAlphabet = (char *)malloc((alphabet_size+1) * sizeof(char));

    memcpy(dAlphabet, alphabet, alphabet_size);
    dAlphabet[alphabet_size] = 0;

    for (i = 0; i < alphabet_size; i++) {
        printf("%c ", alphabet[i]);
        maxRowScore = -INT_MAX;
        for (j = 0; j < pssmLength; j++) {
            pssmScores[alphabet_size * j + i] = scores[pssmLength * i + j];
            printf("%d ", pssmScores[alphabet_size * j + i]);
            nScores++;
        }
        printf("\n");
    }

    for (i = 0; i < pssmLength; i++) {
        maxRowScore = -INT_MAX;
        for (j = 0; j < alphabet_size; j++) {
            if (pssmScores[alphabet_size * i + j] > maxRowScore)
                maxRowScore = pssmScores[alphabet_size * i + j];
        }
        //printf("maxScore: %f maxRowScore: %f\n", maxScore, maxRowScore);
        maxScore += maxRowScore;
    }

    printf("nScores: %d maxScore: %d\n", nScores, maxScore);
    pssm = init_matrix_score(0, pssmLength, alphabet_size, pssmScores, nScores, maxScore / 2);
    pssm->alphabet = dAlphabet;

    free(scores);
    return pssm;
}

void addMinWidthToThresholds(PSSM mat, bwt_width_t *width)
{
    int i;
    for (i = 1; i < mat->length; i++) {
        mat->thresholds[i] += width[i-1].min_drop;
        //fprintf(stderr, "drop[%d] = %f\n", i, width[i].min_drop);
    }
}

void complement_pssm(PSSM mat) {
    int i;
    ubyte_t j;
    int temp_score;
    for (i = 0; i < mat->length; i++) {
        for (j = 0; j < 2; j++) {
            ubyte_t comp_base = 3 - j;
            temp_score = get_score_fast(mat, &j, i);
            mat->scores[mat->offsets[i] + j] = mat->scores[mat->offsets[i] + comp_base];
            mat->scores[mat->offsets[i] + comp_base] = temp_score;
        }
    }
}

void remapAlphabetIndexes(PSSM pssm) {
    int i;
    for (i = 0; i < 256; i++)
        pssm->alphabet_indexes[i] = -1;

    for (i = 0; i < pssm->alphabet_size; i++)  {
        int alph_letter = pssm->alphabet[i];
        pssm->alphabet_indexes[alph_letter] = i;
    }
}


