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

#ifndef PSSM_H_
#define PSSM_H_

#define MAXPSSMSIZE 1024
#define MAX_LINE_LENGTH 1024 
#include "probs.h"

#ifndef BWA_UBYTE
#define BWA_UBYTE
typedef unsigned char ubyte_t;
#endif

struct PSSM
{
int id;
  int length; // can be char
  int order;  // can be char
  int alphabet_size;
  char *alphabet;
  char *alphabet_indexes;
  float *scores;
  int *offsets; //[MAXPSSMSIZE+1];
  float *thresholds; //[MAXPSSMSIZE];
  float *saved_scores; //e[MAXPSSMSIZE+1];
  float *bi; //best inexact match score
  float *be; //best exact match score
  int maxHits;
};

typedef struct PSSM* PSSM;

// Mutators
PSSM init_matrix(int order, int length, int alphabetSize);
PSSM init_matrix_score(int order, int length, int alphabetSize, float *scores, int nScores, float threshold);
void release_matrix(PSSM pssm);
float maximum_score(PSSM pssm, int pos);
void calc_and_set_offsets(PSSM pssm, int order, int length, int alphabetSize);
void calc_and_set_thresholds(PSSM pssm, float threshold);
void calc_and_set_reverse_thresholds(PSSM mat, int start, int end, float threshold);
void set_score(PSSM pssm, const unsigned char *letters, int pos, float score);
static inline void set_length(PSSM pssm, int length) { pssm->length=length; }
void print_horizontalPSSM(PSSM pssm);
void complement_pssm(PSSM mat);

// Acsessors
float get_score(PSSM pssm, const ubyte_t *baseLetter, int pos);
float get_threshold(PSSM pssm, int pos);
float get_global_threshold(PSSM pssm);
int get_length(PSSM pssm);
void calculate_reverse_best_inexact(PSSM mat, float *min_drops);
void calculate_reverse_best_exact(PSSM mat);

// Fast accessors
#define get_threshold_fast(pssm, pos)  (pssm)->thresholds[pos]
#define get_global_threshold_fast(pssm) (pssm)->thresholds[(pssm)->length - 1]
#define get_length_fast(pssm)          (pssm)->length
#define get_num_scores(pssm)           (pssm)->offsets[(pssm)->length]
#define get_score_fast(pssm, base_letter, pos) ( (pssm)->order == 0 ? (pssm)->scores[(pssm)->offsets[pos] + (int)*(base_letter)] : get_score(pssm, base_letter, pos) )
#define get_saved_score(pssm, depth)  (pssm)->saved_scores[depth]
#define set_saved_score(pssm, depth, score)   (pssm)->saved_scores[depth] = score;

#endif
