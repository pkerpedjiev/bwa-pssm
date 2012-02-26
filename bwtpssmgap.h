#ifndef BWTPSSMGAP_H_
#define BWTPSSMGAP_H_

#include "bwt.h"
#include "bwtaln.h"
#include "bwtge.h"
#include "gdsl.h"
//#include "binheap.h"

#ifdef __cplusplus
extern "C" {
#endif

    pssm_heap_t *pssm_init_heap(int max_mm, int max_gapo, int max_gape, const gap_opt_t *opt);
	void pssm_destroy_heap(pssm_heap_t *heap);
	bwt_aln1_t *bwt_match_pssm(bwt_t *const bwt, int len, const ubyte_t *seq, const PSSM mat, bwt_width_t *w, bwt_width_t *seed_w, const gap_opt_t *opt, int *_n_aln, pssm_heap_t *heap);
#ifdef __cplusplus
}
#endif

#endif
