/*
 * seq2pssm.h
 *
 *  Created on: Jun 12, 2009
 *      Author: krogh
 */

#ifndef SEQ2PSSM_H_
#define SEQ2PSSM_H_

#include "shared.h"
#include "probs.h"

/* SEQ2PSSM_H_ */
/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
Probs *alloc_probs(int len, int alphsize);
void free_probs(Probs *P);
void normalize_probs(Probs *P, float pseudocount);
PSSM string_to_pssm(ubyte_t *seq, int len, int alphsize, float match, float mismatch,
			float wcscore);
Probs *alloc_markov_chain(int order, int alphsize);
int sequence_to_pssm(bwa_seq_t *s, int alphsize, float psnp, Probs *mc, float sc_match,
		float sc_mismatch, float sc_wild, int scoretype, float *qualprobs,const gap_opt_t *opt);
float *solexa_ascii_quality_scores();
float *phred_ascii_quality_scores(int base);
float *read_ascii_quality_scores(char *filename);
Probs *qual_to_probs(ubyte_t *seq, ubyte_t *qual, int len, int alphsize,
		float *qualprobs);
void snp_probs(Probs *P, float *q, float psnp);
PSSM prob_to_pssm(Probs *P, Probs *mc);
float highest_scores(PSSM mat);
float mismatch_threshold(PSSM mat, int M);
Probs *markov_chain(bwtint_t *counts, int alphlen);
void set_thresholds(PSSM mat, const gap_opt_t *opt);
/* FUNCTION PROTOTYPES END */

#endif
