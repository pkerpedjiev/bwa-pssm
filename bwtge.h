#include "utils.h"
#include <assert.h>
#include "main.h"
#include "bwt.h"

#include "gdsl/gdsl_interval_heap.h"

#ifndef _BWTGE_H_
#define _BWTGE_H_

struct pssm_entry_st { // recursion stack
	u_int32_t info; // score<<21 | a<<20 | i
	u_int32_t n_mm:8, n_gapo:8, n_gape:8, state:2, n_seed_mm:6;
	bwtint_t k, l; // (k,l) is the SA region of [i,n-1]
    float pssm_score;
    float score_offset; // the difference between the best_possible score and the current score
                        // so pssm_score[i] = mat.be[i] + score_offset
	int last_diff_pos;
    struct pssm_entry_st *next;
    struct pssm_entry_st *prev;
};

typedef struct pssm_entry_st pssm_entry_t;

typedef struct {
    //a pool of pssm_entry_t's that are allocated only once
    pssm_entry_t *entry_list; 

    //the first_empty and last_empty are all drawn from entry_list
    pssm_entry_t *first_empty; 
    pssm_entry_t *last_empty;  
    pssm_entry_t *first_busy;
    pssm_entry_t *last_busy;
    gdsl_interval_heap_t heap;

    //the number of empty entries
    int empty_left;
    int num_busy;

    //the number of available entries
    int size;
} pssm_heap_t;

extern gdsl_element_t 
alloc_pssm_entry (void *pssm_entry);

extern void 
free_pssm_entry (gdsl_element_t e);

extern long int
compare_pssm_entries (gdsl_element_t e1, void* e2);
#endif
