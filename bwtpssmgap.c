#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include "bwtgap.h"
#include "bwtpssmgap.h"
#include "bwtaln.h"
#include "bwtpssm.h"

#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

#define aln_score(m,o,e,p) ((m)*(p)->s_mm + (o)*(p)->s_gapo + (e)*(p)->s_gape)
int debug_print = 0;

float log2f(float arg);

pssm_heap_t *pssm_init_heap(int max_mm, int max_gapo, int max_gape, const gap_opt_t *opt)
{
    pssm_heap_t *gp_heap = (pssm_heap_t *)malloc(sizeof(pssm_heap_t));  
    int i;

    gp_heap->size = opt->max_entries;
    gp_heap->entry_list = (pssm_entry_t *)malloc(gp_heap->size * sizeof(pssm_entry_t));
    gp_heap->num_to_insert = 0;

    gp_heap->heap = gdsl_interval_heap_alloc ("H", alloc_pssm_entry, free_pssm_entry, compare_pssm_entries);
    gdsl_interval_heap_set_max_size(gp_heap->heap, gp_heap->size-15);

    gdsl_interval_heap_flush( gp_heap->heap);

    for (i = 0; i < gp_heap->size - 1; i++) {
        gp_heap->entry_list[i].next = &gp_heap->entry_list[i+1];
    }

    gp_heap->entry_list[gp_heap->size-1].next = 0;

    gp_heap->first_empty = &gp_heap->entry_list[2];
    gp_heap->last_empty = &gp_heap->entry_list[gp_heap->size - 1];
    //gp_heap->empty_left = gp_heap->size;

    gp_heap->first_busy = &gp_heap->entry_list[0];
    gp_heap->last_busy = &gp_heap->entry_list[1];

    gp_heap->first_busy->next = gp_heap->last_busy;

    return gp_heap;
}


void pssm_destroy_heap(pssm_heap_t *gp_heap)
{
    gdsl_interval_heap_free (gp_heap->heap);

    free(gp_heap->entry_list);
    free(gp_heap);
}

/*
static void print_empty(pssm_heap_t *gp_heap)
{
    pssm_entry_t *p = gp_heap->first_empty;
    fprintf(stderr, "e: ");

    while (p != gp_heap->last_empty)  {
        fprintf(stderr, "%p ", p);
        p = p->next;
    }
    fprintf(stderr, "%p ", p);

    fprintf(stderr, "\n");
}

static void check_empty_integrity(pssm_heap_t *gp_heap)
{
    pssm_entry_t *p;

    int count = 1;

    //fprintf(stderr, "update_heap: gp_heap->first_empty: %x gp_heap->last_empty: %x\n", gp_heap->first_empty, gp_heap->last_empty);

    p = gp_heap->first_empty;
    while (p != gp_heap->last_empty) {
        //fprintf(stderr, "%x ", p);
        p  = p->next;
        count++;
    }
    //fprintf(stderr, "%x ", p);
    //fprintf(stderr, "\n");

    fprintf(stderr, "count: %d, empty left: %d\n", count, gp_heap->empty_left);
    //print_empty(gp_heap);
    //
    assert(count == gp_heap->empty_left);
    assert(gp_heap->first_empty != gp_heap->last_empty);
}

*/
static void remove_from_empty_list(pssm_heap_t *gp_heap, pssm_entry_t *p)
{
    //fprintf(stderr, "returning to empty list %x\n", p);

    p->prev->next = p->next;
    p->next->prev = p->prev;

    gp_heap->last_empty->next = p;
    gp_heap->last_empty = p;
    gp_heap->empty_left++;

    //check_empty_integrity(gp_heap);
}

static inline pssm_entry_t *update_heap(pssm_heap_t *gp_heap)
{

    pssm_entry_t *last_busy = gp_heap->last_busy;
    pssm_entry_t *p;


    // take an entry from the empty list
    p = gp_heap->first_empty;
    gp_heap->first_empty = p->next;
    gp_heap->empty_left--;

    // insert that entry into the busy list
    last_busy->prev->next = p;
    p->prev = last_busy->prev;
    p->next = last_busy;
    last_busy->prev = p;

        //fprintf(stderr, "taking from empty list: %x\n", p);
    //check_empty_integrity(gp_heap);

    return p;
}

static void gap_push(pssm_heap_t *gp_heap, int id, int a, int i, bwtint_t k, bwtint_t l, int n_mm, int n_gapo, int n_gape,
        int state, int is_diff, const gap_opt_t *opt, int pssm_score, int score_offset)
{
    int score;
    pssm_entry_t *p;
    id=0;

    assert(gp_heap->empty_left > 0);
    p = update_heap(gp_heap);

    score = aln_score(n_mm, n_gapo, n_gape, opt);
    p->info = (u_int32_t)score<<21 | a<<20 | i; p->k = k; p->l = l;
    p->n_mm = n_mm; p->n_gapo = n_gapo; p->n_gape = n_gape; p->state = state;
    p->last_diff_pos=0;
    p->score_offset = score_offset;
    p->pssm_score = pssm_score;
    if (is_diff) p->last_diff_pos = i;
    //fprintf(stderr, "inserting score: %d\n", score_offset);
    //fprintf(stderr, "score: %f &p->score_offset: %x (p): %x p+n: %x\n", *(int *)(((void *)p) + 28), &p->score_offset, p, ((void *)p) + 28);
    //fprintf(stderr, "a - b: %d\n", (int)&p->score_offset - (int)p);
    //

    gp_heap->to_insert[gp_heap->num_to_insert++] = p;
    /*
    p1 = gdsl_interval_heap_insert(gp_heap->heap, (void *)p);
    //gdsl_raw_heap_dump(p);
    if (p1 != NULL)
        remove_from_empty_list(gp_heap, p1);
        */
}

static void gap_finish_push(pssm_heap_t *gp_heap) {
    int i;
    pssm_entry_t *p1;

    //fprintf(stderr, "finish pushing: %d\n", gp_heap->num_to_insert);
    for (i = 0; i < gp_heap->num_to_insert; i++) {
        p1 = NULL;
        //fprintf(stderr, "inserting: %x score: %d\n", gp_heap->to_insert[i], gp_heap->to_insert[i]->score_offset);
        p1 = gdsl_interval_heap_insert(gp_heap->heap, gp_heap->to_insert[i]);
        if (p1 != NULL) {
            //fprintf(stderr, "removing %x score: %d\n", p1, p1->score_offset);
            remove_from_empty_list(gp_heap, p1);
        }
    }

    gp_heap->num_to_insert=0;
}

static void gap_pop(pssm_heap_t *gp_heap, int id, pssm_entry_t *e)
{
    pssm_entry_t *p = (pssm_entry_t *)gdsl_interval_heap_remove_max (gp_heap->heap);

    memcpy(e, p, sizeof(pssm_entry_t));
    //assume p is part of the busy list and is between the first and last
    //
    remove_from_empty_list(gp_heap, p);
}

/*
static int gap_destroy_min(pssm_heap_t *gp_heap) {

    pssm_entry_t *p = (pssm_entry_t *)gdsl_interval_heap_remove_min (gp_heap->heap);
    //print_empty(gp_heap);
    //fprintf(stderr, "destroy_min: gp_heap->first_empty: %x gp_heap->last_empty: %x p: %x\n", gp_heap->first_empty, gp_heap->last_empty, p);
    remove_from_empty_list(gp_heap, p);
    //print_empty(gp_heap);
    //fprintf(stderr, "end destroy_min\n");
    return p->score_offset;
}
*/

static void gap_reset_heap(pssm_heap_t *gp_heap)
{
    gdsl_interval_heap_flush( gp_heap->heap);
    //gp_heap->heap = gdsl_interval_heap_alloc ("H", alloc_pssm_entry, free_pssm_entry, compare_pssm_entries);
    //gdsl_interval_heap_set_max_size(gp_heap->heap, gp_heap->size-10);

    if (gp_heap->first_busy->next != gp_heap->last_busy) {
        gp_heap->last_empty->next = gp_heap->first_busy->next;
        gp_heap->last_empty = gp_heap->last_busy->prev;
    }

    gp_heap->first_busy->next = gp_heap->last_busy;
    gp_heap->last_busy->prev = gp_heap->first_busy;
    gp_heap->empty_left = gp_heap->size-2;
}


static inline void gap_shadow(int x, int len, bwtint_t max, int last_diff_pos, bwt_width_t *w)
{
    int i, j;
    for (i = j = 0; i < last_diff_pos; ++i) {
        if (w[i].w > x) w[i].w -= x;
        else if (w[i].w == x) {
            w[i].bid = 1;
            w[i].w = max - (++j);
        } // else should not happen
    }
}

static inline int int_log2(uint32_t v)
{
    int c = 0;
    if (v & 0xffff0000u) { v >>= 16; c |= 16; }
    if (v & 0xff00) { v >>= 8; c |= 8; }
    if (v & 0xf0) { v >>= 4; c |= 4; }
    if (v & 0xc) { v >>= 2; c |= 2; }
    if (v & 0x2) c |= 1;
    return c;
}

bwt_aln1_t *bwt_match_pssm(bwt_t *const bwt, int len, const ubyte_t *seq, const PSSM mat, bwt_width_t *width, bwt_width_t *seed_width, const gap_opt_t *opt, int *_n_aln, pssm_heap_t *gp_heap)
{
    int best_score = aln_score(opt->max_diff+1, opt->max_gapo+1, opt->max_gape+1, opt);
    int best_diff = opt->max_diff + 1, max_diff = opt->max_diff;
    int best_cnt = 0;
    int j, _j, n_aln, m_aln;
    int visited = 0;
    int gap_open_penalty = -(int)(1000 *opt->p_gapo); 
    int gap_extension_penalty = -(int)(1000 * opt->p_gape);
    int deletion_penalty = -(int)(1000 * opt->p_del);
    bwt_aln1_t *aln;
    int hit_found = 0;
    int num_hits_found = 0;
    int best_found = -INT_MAX;
    int min_score = -INT_MAX;
    int desired_mapq = 5000;
    int curr_offset;

    m_aln = 4; n_aln = 0;
    aln = (bwt_aln1_t*)calloc(m_aln, sizeof(bwt_aln1_t));

    // check whether there are too many N
    for (j = _j = 0; j < len; ++j)
        if (seq[j] > 3) ++_j;
    if (_j > max_diff) {
        *_n_aln = n_aln;
        return aln;
    }

    gap_reset_heap(gp_heap); // reset heap
    gap_push(gp_heap, mat->id, 0, len, 0, bwt->seq_len, 0, 0, 0, 0, 0, opt, 0.0, 0.0);
    gap_finish_push(gp_heap);

    //the most entries that can be pushed onto the heap is at most 7:
    //one match, three mismatches, open, extension, deletion
    while (gp_heap->empty_left < (gp_heap->size-2) ) {
        pssm_entry_t e;
        int a, i, m, m_seed = 0, allow_diff, allow_M, tmp;
        bwtint_t k, l, cnt_k[4], cnt_l[4], occ;
        int curr_score;
        int curr_threshold;
        min_score = -INT_MAX;

        g_visited++;
        visited++;

        gap_pop(gp_heap, mat->id, &e); // get the best entry

        k = e.k; l = e.l; // SA interval
        a = e.info>>20&1; i = e.info&0xffff; // strand, length

        if (!(opt->mode & BWA_MODE_NONSTOP) && best_found > mat->be[mat->length-1] + e.score_offset + desired_mapq) {
            break;
        }

        int max_entries = 0;
        //fprintf(stderr, "pssm #1 id:%d %d \t[%d][%d,%d,%d,%d,%c]\t[%d,%d,%d]\t[%u,%lu]\t[%lu,%lu]\t%d\t[%6d, **%6d**, %6d, %6d]\n", mat->id, i, max_entries, gp_heap->empty_left, a, i, seq[i], "MID"[e.state], e.n_mm,     e.n_gapo, e.n_gape, width[i-1].bid, width[i-1].w, k, l, e.last_diff_pos, curr_score, e.score_offset, mat->thresholds[i], mat->bi[i]);

        m = max_diff - (e.n_mm + e.n_gapo);
        if (i == mat->length)
            curr_score = 0;
        else 
            curr_score = mat->be[mat->length-1] - mat->be[i] + e.score_offset;

        if (opt->mode & BWA_MODE_GAPE) m -= e.n_gape;
        if (m < 0) {
            continue;
        }
        if (seed_width) { // apply seeding
            m_seed = opt->max_seed_diff - (e.n_mm + e.n_gapo);
            if (opt->mode & BWA_MODE_GAPE) m_seed -= e.n_gape;
        }

        // check whether a hit is found
        hit_found = 0;
        if (i == 0) {
            hit_found = 1;
            num_hits_found += 1;
            if (i == 0) {
                if (curr_score > best_found) {
                    if (num_hits_found >= 2) {
                        calc_and_set_reverse_thresholds(mat, 1, get_length(mat), curr_score);
                        addMinWidthToThresholds(mat, width);
                    }
                    best_found = curr_score;
                }
            }
        }
        else if (curr_score + mat->bi[i-1] < mat->thresholds[0]) {
            if (bwt_match_exact_alt(bwt, i, seq, &k, &l)) {
                hit_found = 1;
                num_hits_found += 1;
                e.score_offset = 0;
            }
            else {
                continue; // no hit, skip
            }
        }

        if (hit_found) { // action for found hits
            int score = aln_score(e.n_mm, e.n_gapo, e.n_gape, opt);

            int do_add = 1;
            if (n_aln == 0) {
                best_score = score;
                best_diff = e.n_mm + e.n_gapo;
                if (opt->mode & BWA_MODE_GAPE) best_diff += e.n_gape;
                if (!(opt->mode & BWA_MODE_NONSTOP))
                    max_diff = (best_diff + 1 > opt->max_diff)? opt->max_diff : best_diff + 1; // top2 behaviour
            }
            if (score == best_score) best_cnt += l - k + 1;
            else if (best_cnt > opt->max_top2) {
                break; // top2b behaviour
            }
            if (e.n_gapo) { // check whether the hit has been found. this may happen when a gap occurs in a tandem repeat
                for (j = 0; j != n_aln; ++j)
                    if (aln[j].k == k && aln[j].l == l) {
                        break;
                    }
                if (j < n_aln) do_add = 0;
            }
            if (do_add) { // append
                bwt_aln1_t *p;
                gap_shadow(l - k + 1, len, bwt->seq_len, e.last_diff_pos, width);
                if (n_aln == m_aln) {
                    m_aln <<= 1;
                    aln = (bwt_aln1_t*)realloc(aln, m_aln * sizeof(bwt_aln1_t));
                    memset(aln + m_aln/2, 0, m_aln/2*sizeof(bwt_aln1_t));
                }
                p = aln + n_aln;
                p->n_mm = e.n_mm; p->n_gapo = e.n_gapo; p->n_gape = e.n_gape; p->a = a;
                p->k = k; p->l = l;
                p->score = score;
                p->pssm_score = curr_score;
                p->pssm = 1;
                ++n_aln;
            }
            continue;
        }

        --i;
        bwt_2occ4(bwt, k - 1, l, cnt_k, cnt_l); // retrieve Occ values
        occ = l - k + 1;
        // test whether diff is allowed
        allow_diff = allow_M = 1;
        curr_threshold = mat->thresholds[i];

        if (i > 0) {
            int ii = i - (len - opt->seed_len);
            if (width[i-1].bid > m-1) allow_diff = 0;
            else if (width[i-1].bid == m-1 && width[i].bid == m-1 && width[i-1].w == width[i].w) allow_M = 0;
            if (seed_width && ii > 0) {
                if (seed_width[ii-1].bid > m_seed-1) allow_diff = 0;
                else if (seed_width[ii-1].bid == m_seed-1 && seed_width[ii].bid == m_seed-1
                        && seed_width[ii-1].w == seed_width[ii].w) allow_M = 0;
            }
        }
        // indels
        tmp = (opt->mode & BWA_MODE_LOGGAP)? int_log2(e.n_gape + e.n_gapo)/2+1 : e.n_gapo + e.n_gape;
        if (allow_diff && i >= opt->indel_end_skip + tmp && len - i >= opt->indel_end_skip + tmp) {
            switch (e.state) {
                case STATE_M: // gap open
                    if (e.n_gapo < opt->max_gapo) { // gap open is allowed
                        curr_offset = e.score_offset - gap_open_penalty;

                        // insertion
                        if ((curr_score - gap_open_penalty > curr_threshold) & 
                            (curr_offset > min_score))  {
                                gap_push(gp_heap, mat->id, a, i, k, l, e.n_mm, e.n_gapo + 1, e.n_gape, STATE_I, 1, opt, curr_score - gap_open_penalty, curr_offset);
                        }

                        curr_offset = e.score_offset - deletion_penalty;
                            // deletion
                        if (curr_score - deletion_penalty > curr_threshold && curr_offset > min_score) { 
                            for (j = 0; j != 4; ++j) {
                                k = bwt->L2[j] + cnt_k[j] + 1;
                                l = bwt->L2[j] + cnt_l[j];
                                if (k <= l) 
                                    gap_push(gp_heap, mat->id, a, i + 1, k, l, e.n_mm, e.n_gapo + 1, e.n_gape, STATE_D, 1, opt, curr_score - deletion_penalty, curr_offset);
                            }
                        }
                    }
                    break;
                case STATE_I: //extension of an insertion
                    if (e.n_gape < opt->max_gape) { // gap extention is allowed
                            curr_offset = e.score_offset - gap_extension_penalty;
                            if (curr_score - gap_extension_penalty > curr_threshold && curr_offset > min_score) {
                                gap_push(gp_heap, mat->id, a, i, k, l, e.n_mm, e.n_gapo, e.n_gape + 1, STATE_I, 1, opt, curr_score - gap_extension_penalty, curr_offset);
                        }
                    }
                    break;
                case STATE_D: //extension of a deletion
                    if (e.n_gape < opt->max_gape) { // gap extention is allowed
                        if (e.n_gape + e.n_gapo < max_diff || occ < opt->max_del_occ) {
                            curr_offset = e.score_offset - gap_extension_penalty;
                            if (curr_score - gap_extension_penalty > curr_threshold && curr_offset > min_score) {
                                for (j = 0; j != 4; ++j) {
                                    k = bwt->L2[j] + cnt_k[j] + 1;
                                    l = bwt->L2[j] + cnt_l[j];
                                    
                                    if (k <= l) gap_push(gp_heap, mat->id, a, i + 1, k, l, e.n_mm, e.n_gapo, e.n_gape + 1, STATE_D, 1, opt, curr_score - gap_extension_penalty, curr_offset);
                                }
                            }
                        }
                    }
                break;
                }
            }

        // mismatches
        if (allow_diff && allow_M) { // mismatch is allowed
            int some_score;
            if (i == 0) {
                some_score = mat->be[0] + e.score_offset;
            } else {
                some_score = -mat->be[i] + mat->be[i-1] + e.score_offset;
            }
            for (j = 4; j >= 1; --j) {
                ubyte_t c = (seq[i] + j) & 3;
                int is_mm = (j != 4 || seq[i] > 3);
                int base_score = get_score_fast(mat, &c, i);
                int curr_offset = some_score + base_score;
                //int score = -((mat->be[i] - mat->be[i-1]) - base_score) + e.score_offset;
                
                if (curr_score + base_score < curr_threshold && curr_offset < min_score)  {
                    continue;
                }

                k = bwt->L2[c] + cnt_k[c] + 1;
                l = bwt->L2[c] + cnt_l[c];
                if (k <= l) {
                    gap_push(gp_heap, mat->id, a, i, k, l, e.n_mm + is_mm, e.n_gapo, e.n_gape, STATE_M, is_mm, opt, curr_score + base_score, curr_offset);
                }
            }
        } else if (seq[i] < 4) { // try exact match only
            ubyte_t c = seq[i] & 3;
            int base_score = get_score_fast(mat, &c, i);
            int curr_offset;
            if (i == 0) {
                curr_offset = -((mat->be[0]) - base_score) + e.score_offset;
            } else {
                curr_offset = -((mat->be[i] - mat->be[i-1]) - base_score) + e.score_offset;
            }

            if (curr_offset > min_score) {

            k = bwt->L2[c] + cnt_k[c] + 1;
            l = bwt->L2[c] + cnt_l[c];


                if (k <= l) gap_push(gp_heap, mat->id, a, i, k, l, e.n_mm, e.n_gapo, e.n_gape, STATE_M, 0, opt, curr_score + base_score, curr_offset);
            }
        }
        gap_finish_push(gp_heap);
    }

    *_n_aln = n_aln;
    return aln;
    }
