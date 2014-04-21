#ifndef BWTPSSM_H
#define BWTPSSM_H

#include <stdint.h>
#include "pssm.h"
#include "bwt.h"

#ifdef __cplusplus
extern "C" {
#endif

	void bwa_pssm_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt);

    bwa_seq_t *bwa_read_pssm_seq(bwa_seqio_t *bs, int n_needed, int *n, int mode, int trim_qual, Probs *mc, float *qualprobs,const gap_opt_t *opt);
	void bwa_cal_pssm_sa_reg_gap(int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt);
    bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);

	/* rgoya: Temporary clone of aln_path2cigar to accomodate for bwa_cigar_t,
	__cigar_op and __cigar_len while keeping stdaln stand alone */
#include "stdaln.h"

    void addMinWidthToThresholds(PSSM mat, bwt_width_t *width);
    Probs *markov_chain(bwtint_t *counts, int alphlen);
    float *phred_ascii_quality_scores(int base);
#ifdef __cplusplus
}
#endif

#endif
