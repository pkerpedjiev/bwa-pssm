#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "pssm.h"
#include "bwtaln.h"
#include "bwtgap.h"
#include "bwtpssm.h"
#include "bwtpssmgap.h"
#include "utils.h"
#include "seq2pssm.h"

#ifdef HAVE_PTHREAD
#define THREAD_BLOCK_SIZE 1024
#include <pthread.h>
#endif

float log2f(float arg);

/* Calculate the difference between the greatest score and the least score in at a
   location in the pssm */
double calc_min_drop_at_pos(PSSM mat, int pos) {
    double score;
    double best_score = -DBL_MAX;
    double second_best_score = -DBL_MAX;
    int i=pos;
    ubyte_t j;

    //four letters in the alphabet

    for (j = 0; j < 4; j++) {
        score = get_score_fast(mat, &j, i);
        //printf(stderr, "score: %f i: %d j: %d\n", score, i, j);
        if (score > best_score) {
            best_score = score;
        }
    }

    for (j = 0; j < 4; j++) {
        score = get_score_fast(mat, &j, i);
        if (score < best_score && score > second_best_score) {
            second_best_score = score;
        }
    }


    return best_score - second_best_score;
}

/* Calculate the minimum score drop in the first n characters after start of the pssm.
 * This that is calculate the difference between best_score[i] - 2nd_best_score[i]
 * for each start <= i <= n
 */
float calc_min_drop(int *min_drops, int start, int end) {
    int i;
    double min_drop = DBL_MAX;

    /*
       if (n > end - start) {
    //uh-oh, this should never be
    return DBL_MAX;
    }
     */

    for (i = start; i <= end; i++) 
        if (min_drops[i] < min_drop) 
            min_drop = min_drops[i];

    return min_drop;
}


// width must be filled as zero
static int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, int *min_drops, bwt_width_t *width)
{   
    bwtint_t k, l, ok, ol;
    int i, bid, start = 0;
    float min_drop = 0;

    bid = 0;
    k = 0; l = bwt->seq_len;
    for (i = 0; i < len; ++i) {
        ubyte_t c = str[i];
        if (c < 4) {
            bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
            k = bwt->L2[c] + ok + 1;
            l = bwt->L2[c] + ol;
        }
        if (k > l || c > 3) { // then restart
            k = 0;
            l = bwt->seq_len;
            ++bid;
            min_drop += calc_min_drop(min_drops, start, i);
            start = i;
        }
        width[i].w = l - k + 1;
        width[i].bid = bid;
        width[i].min_drop = min_drop;
    }
    width[len].w = 0;
    width[len].bid = ++bid;
    width[len].min_drop = min_drop;
    return bid;
}

float min(float a, float b) {
    return a > b ? b : a;
}

void bwa_cal_pssm_sa_reg_gap(int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt)
{
    int i, j, max_l = 0, max_len;
    bwt_width_t *w, *seed_w;
    int *min_drops;
    gap_opt_t local_opt = *opt;
    pssm_heap_t *gp_heap;

    // initiate priority stack
    for (i = max_len = 0; i != n_seqs; ++i)
        if (seqs[i].len > max_len) max_len = seqs[i].len;
    if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
    if (local_opt.max_diff < local_opt.max_gapo) local_opt.max_gapo = local_opt.max_diff;

    if (!local_opt.pssm_ratio_provided) {
        local_opt.pssm_ratio = ((float)local_opt.max_diff) - local_opt.pssm_ratio_discount;
        //local_opt.max_diff += 1;
    }

    local_opt.max_diff = 30;

    gp_heap = pssm_init_heap(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);

    seed_w = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));

    w = 0;
    for (i = 0; i != n_seqs; ++i) {
        bwa_seq_t *p = seqs + i;
#ifdef HAVE_PTHREAD
        if (i % opt->n_threads != tid) continue;
#endif
        //fprintf(stderr, "i: %d tid: %d\n", i, tid);
        p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
        min_drops = (float *)calloc(p->len+1, sizeof(float));
        p->mat->id = i;

        for (j = 0; j < p->len; j++) {
                float curr_min_drop = calc_min_drop_at_pos(p->mat, j);
                
                curr_min_drop = min(curr_min_drop, -opt->p_gapo);
                curr_min_drop = min(curr_min_drop, -opt->p_gape);
                curr_min_drop = min(curr_min_drop, -opt->p_del);
               
                min_drops[j] = curr_min_drop;
            }

        if (max_l < p->len) {
            max_l = p->len;
            w = (bwt_width_t*)realloc(w, (max_l + 1) * sizeof(bwt_width_t));
            memset(w, 0, (max_l + 1) * sizeof(bwt_width_t));
        }

        bwt_cal_width(bwt, p->len, p->seq, min_drops, w);
        if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
        if (!opt->pssm_ratio_provided) {
            local_opt.pssm_ratio = ((float)local_opt.max_diff) - opt->pssm_ratio_discount;
            //local_opt.max_diff += 1;
        }
        local_opt.max_diff = 30;

        //fprintf(stderr, "local_opt.pssm_ratio: %f\n", local_opt.pssm_ratio);
        local_opt.seed_len = opt->seed_len < p->len? opt->seed_len : 0x7fffffff;
        if (p->len > opt->seed_len) {
            bwt_cal_width(bwt, opt->seed_len, p->seq + (p->len - opt->seed_len), min_drops, seed_w);
        }

        set_thresholds(p->mat, &local_opt);
        calculate_reverse_best_inexact(p->mat, min_drops);
        calculate_reverse_best_exact(p->mat);
        //local_opt.max_diff=30;
        addMinWidthToThresholds(p->mat, w);
        complement_pssm(p->mat);

        //core function
        for (j = 0; j < p->len; ++j) // we need to complement
            p->seq[j] = p->seq[j] > 3? 4 : 3 - p->seq[j];

            /*
        for (j = 0; j < p->len; ++j) {
            float max_score = -100.0;
            for (k = 0; k < 4; k++) {
                if (get_score_fast(p->mat, &k, j) > max_score) {
                    p->seq[j] = k;
                    max_score = get_score_fast(p->mat, &k, j);
                }
            }
        }
        */

        //fprintf(stderr, "max_diff: %d pssm_ratio: %f\n", local_opt.max_diff, local_opt.pssm_ratio);
        p->aln = bwt_match_pssm(bwt, p->len, p->seq, p->mat, w, p->len <= opt->seed_len? 0 : seed_w, &local_opt, &p->n_aln, gp_heap);
        // store the alignment
        free(min_drops); 
        free(p->name); free(p->seq); free(p->rseq); free(p->qual); 
        free(p->rqual);
        release_matrix(p->mat);
        p->name = 0; p->seq = p->rseq = p->qual = p->rqual = 0;
    }
    free(seed_w); free(w); 
    pssm_destroy_heap(gp_heap);
}

#ifdef HAVE_PTHREAD
typedef struct {
    int tid;
    bwt_t *bwt;
    int n_seqs;
    bwa_seq_t *seqs;
    const gap_opt_t *opt;
} thread_aux_t;

static void *pssm_worker(void *data)
{
    thread_aux_t *d = (thread_aux_t*)data;
    bwa_cal_pssm_sa_reg_gap(d->tid, d->bwt, d->n_seqs, d->seqs, d->opt);
    return 0;
}
#endif

void bwa_pssm_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt)
{
    int i, n_seqs, tot_seqs = 0;
    bwa_seq_t *seqs;
    bwa_seqio_t *ks;
    clock_t t;
    bwt_t *bwt;
    float *qualprobs;
    int qbase=33;
    Probs *mc;

    //calculate the quality score probabilities
    qualprobs = phred_ascii_quality_scores(qbase);

    // initialization
    g_visited = 0;
    ks = bwa_open_reads(opt->mode, fn_fa);

    { // load BWT
        char *str = (char*)calloc(strlen(prefix) + 10, 1);
        //strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
        strcpy(str, prefix); strcat(str, ".bwt"); bwt = bwt_restore_bwt(str);
        free(str);
    }

    mc = markov_chain(bwt->L2, 4);

    // core loop
    fwrite(opt, sizeof(gap_opt_t), 1, stdout);
    while ((seqs = bwa_read_pssm_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual, mc, qualprobs, opt)) != 0) {
        tot_seqs += n_seqs;
        t = clock();


#ifdef HAVE_PTHREAD
        if (opt->n_threads <= 1) { // no multi-threading at all
            bwa_cal_pssm_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
        } else {
            pthread_t *tid;
            pthread_attr_t attr;
            thread_aux_t *data;
            int j;
            pthread_attr_init(&attr);
            pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
            tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
            for (j = 0; j < opt->n_threads; ++j) {
                data[j].tid = j; data[j].bwt = bwt; 
                data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
                pthread_create(&tid[j], &attr, pssm_worker, data + j);
            }
            for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
            free(data); free(tid);
        }
#else
        bwa_cal_pssm_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
#endif

        fprintf(stderr, "[bwa_pssm_core] calculate SA coordinate... ");
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

        t = clock();
        fprintf(stderr, "[bwa_pssm_core] write to the disk... ");
        for (i = 0; i < n_seqs; ++i) {
            bwa_seq_t *p = seqs + i;
            fwrite(&p->n_aln, 4, 1, stdout);
            if (p->n_aln) fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, stdout);
        }
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

        bwa_free_read_seq(n_seqs, seqs);
        fprintf(stderr, "[bwa_pssm_core] %d sequences have been processed.\n", tot_seqs);
    }
    fprintf(stderr, "g_visited: %lu\n", g_visited);

    // destroy
    bwt_destroy(bwt);
    //bwt_destroy(bwt[0]);
    bwa_seq_close(ks);
}

int bwa_pssm(int argc, char *argv[])
{
    int c, opte = -1;
    gap_opt_t *opt;

    opt = gap_init_opt();

    opt->max_entries = 400;
    
    while ((c = getopt(argc, argv, "pn:z:y:o:e:i:d:l:k:cLR:m:t:NM:O:E:D:G:P:q:f:b012IB:")) >= 0) {
        switch (c) {
            case 'n':
                if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
                else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
                break;
            case 'p': opt->parclip = 1; break;
            case 'z': opt->pssm_ratio = atof(optarg); opt->pssm_ratio_provided=1; break;
            case 'y': opt->pssm_ratio_discount = atof(optarg); break;
            case 'o': opt->max_gapo = atoi(optarg); break;
                      //case 'D': opt->debug = 1; break;
            case 'e': opte = atoi(optarg); break;
            case 'M': opt->s_mm = atoi(optarg); break;
            case 'O': opt->p_gapo = log2f(atof(optarg)); break;
            case 'E': opt->p_gape = log2f(atof(optarg)); break;
            case 'D': opt->p_del = log2f(atof(optarg)); break;
            case 'G': load_error_model(opt->error_lookup, optarg); opt->use_error_model=1; break;
            case 'P': opt->prior = atof(optarg); break;
            case 'd': opt->max_del_occ = atoi(optarg); break;
            case 'i': opt->indel_end_skip = atoi(optarg); break;
            case 'l': opt->seed_len = atoi(optarg); break;
            case 'k': opt->max_seed_diff = atoi(optarg); break;
            case 'm': opt->max_entries = atoi(optarg); break;
            case 't': opt->n_threads = atoi(optarg); break;
            case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
            case 'R': opt->max_top2 = atoi(optarg); break;
            case 'q': opt->trim_qual = atoi(optarg); break;
            case 'c': opt->mode &= ~BWA_MODE_COMPREAD; break;
            case 'N': opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff; break;
            case 'f': xreopen(optarg, "wb", stdout); break;
            case 'b': opt->mode |= BWA_MODE_BAM; break;
            case '0': opt->mode |= BWA_MODE_BAM_SE; break;
            case '1': opt->mode |= BWA_MODE_BAM_READ1; break;
            case '2': opt->mode |= BWA_MODE_BAM_READ2; break;
            case 'I': opt->mode |= BWA_MODE_IL13; break;
            case 'Y': opt->mode |= BWA_MODE_CFY; break;
            case 'B': opt->mode |= atoi(optarg) << 24; break;
            default: return 1;
        }
    }
    if (opte > 0) {
        opt->max_gape = opte;
        opt->mode &= ~BWA_MODE_GAPE;
    }


    if (optind + 2 > argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   bwa pssm [options] <prefix> <in.fq>\n\n");
        fprintf(stderr, "Options: -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n",
                BWA_AVG_ERR, opt->fnr);
        fprintf(stderr, "         -o INT    maximum number (int) or fraction of gap opens [%d]\n", opt->max_gapo);
        fprintf(stderr, "         -e INT    maximum number (int) or fraction of gap extensions, -1 for disabling long gaps [-1]\n");
        fprintf(stderr, "         -i INT    do not put an indel within INT bp towards the ends [%d]\n", opt->indel_end_skip);
        fprintf(stderr, "         -d INT    maximum occurrences for extending a long deletion [%d]\n", opt->max_del_occ);
        fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
        fprintf(stderr, "         -k INT    maximum differences in the seed [%d]\n", opt->max_seed_diff);
        fprintf(stderr, "         -m INT    maximum entries in the queue [%d]\n", opt->max_entries);
        fprintf(stderr, "         -t INT    number of threads [%d]\n", opt->n_threads);
        fprintf(stderr, "         -M INT    mismatch penalty [%d]\n", opt->s_mm);
        fprintf(stderr, "         -O INT    gap open penalty [%d]\n", opt->s_gapo);
        fprintf(stderr, "         -E INT    gap extension penalty [%d]\n", opt->s_gape);
        fprintf(stderr, "         -P INT    posterior probability of an alignment [%f]\n", opt->prior);
        fprintf(stderr, "         -G INT    error model table\n");
        fprintf(stderr, "         -R INT    stop searching when there are >INT equally best hits [%d]\n", opt->max_top2);
        fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
        fprintf(stderr, "         -f FILE   file to write output to instead of stdout\n");
        fprintf(stderr, "         -B INT    length of barcode\n");
        fprintf(stderr, "         -c        input sequences are in the color space\n");
        fprintf(stderr, "         -L        log-scaled gap penalty for long deletions\n");
        fprintf(stderr, "         -N        non-iterative mode: search for all n-difference hits (slooow)\n");
        fprintf(stderr, "         -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
        fprintf(stderr, "         -b        the input read file is in the BAM format\n");
        fprintf(stderr, "         -0        use single-end reads only (effective with -b)\n");
        fprintf(stderr, "         -1        use the 1st read in a pair (effective with -b)\n");
        fprintf(stderr, "         -2        use the 2nd read in a pair (effective with -b)\n");
        fprintf(stderr, "\n");
        return 1;
    }

    if (opt->prior > 1.0 || opt->prior <= 0) {
        fprintf(stderr, "Invalid prior value of %f. Using 0.8 instead.\n", opt->prior);
        opt->prior = 0.8;
    }

    if (opt->fnr > 0.0) {
        int i, k;
        for (i = 17, k = 0; i <= 250; ++i) {
            int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, opt->fnr);
            if (l != k) fprintf(stderr, "[bwa_aln] %dbp reads: max_diff = %d\n", i, l);
            k = l;
        }
    }
    bwa_pssm_core(argv[optind], argv[optind+1], opt);
    free(opt);
    return 0;
}


