/*

 * seq2pssm.c
 *
 *  Created on: Jun 12, 2009
 *      Author: krogh
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <limits.h>

#include "bwtaln.h"
#include "pssm.h"
#include "probs.h"
#include "seq2pssm.h"
//#include "bitsandpieces.h"

unsigned char nst_nt4_table[256];

float log2f(float arg);
// #define DEBUG

/* Does not allocate the count array and the order array, since they
   are only used in special occasions
*/
Probs *alloc_probs(int len, int alphsize) {
  Probs *P = (Probs *)malloc(sizeof(Probs));
  P->len = len;
  P->alphsize = alphsize;
  P->order = 0;
  P->powers = NULL;
  P->counts = NULL;
  P->p = NULL;
  if (len) P->p = (float *)calloc((len*(alphsize+1)),sizeof(float));
  return P;
}
void free_probs(Probs *P) {
  if (P->powers) free(P->powers);
  if (P->counts) free(P->counts);
  if (P->p) free(P->p);
  free(P);
}



#ifdef DEBUG
static void print_probs(FILE *fp, Probs *P) {
	int i, b;
	/* If order is 0, it is assumed that it is len x alphsize+1 probabilities */
	if (P->order==0) {
		for (i=0; i<P->len; ++i) fprintf(fp,"  %8d",i+1);
		fprintf(fp,"\nPROBS ---------\n");
		if (P->p) {
			for (b=0; b<=P->alphsize; ++b) {
				for (i=b; i<P->len*(P->alphsize+1); i+=P->alphsize+1) fprintf(fp,"  %8.6f",P->p[i]);
				fprintf(fp,"\n");
			}
		}
		fprintf(fp,"\nCOUNTS ---------\n");
		if (P->counts) {
			for (b=0; b<=P->alphsize; ++b) {
				for (i=b; i<P->len*(P->alphsize+1); i+=P->alphsize+1) fprintf(fp,"  %8d",P->counts[i]);
				fprintf(fp,"\n");
			}
		}
	}
	/* Assume a Markov Chain */
	else {
		int k, lett[256];
		for (k=0;k<=P->order;++k) lett[k]=0;

		for (i=0; i< P->len*(P->alphsize+1); ++i) {
			if ( i%(P->alphsize+1)==0 && i>0) fprintf(fp,"\n");
			else if (i>0) fprintf(fp,"  ");
			for (k=0;k<=P->order;++k) fprintf(fp,"%1d",lett[k]);
			fprintf(fp,":%8.6f",P->p[i]);
			if (P->counts) fprintf(fp,":%8d",P->counts[i]);
			/* Keep track of the word */
			for (k=P->order; k>=0; --k) {
				lett[k] +=1;
				if (lett[k] > P->alphsize) lett[k]=0;
				else break;
			}
		}
		fprintf(fp,"\n");
	}
}
#endif



/* Make a normalized probability vector (P->p)
 * If counts are given use these, otherwise use P->p
 * Add psudocount before normalization
 * Prob of wildcard = 1/alphsize
*/
void normalize_probs(Probs *P, float pseudocount) {
	int n, i, L;
	float sum, *p;

	p=P->p;
	L = P->len*(P->alphsize+1);
	if (P->counts) for (n=0; n<L;++n) p[n]=P->counts[n];
	for (n=0; n<L;++n) p[n] += pseudocount;
	for (n=0; n<P->len;++n) {
		sum=0.;
		for (i=0; i<P->alphsize; ++i) sum+=p[i];
		if (sum>1.e-50) for (i=0; i<P->alphsize; ++i) p[i] = p[i]/sum;
		p[i] = 1.0/P->alphsize;
		p+=P->alphsize+1;
	}
#ifdef DEBUG
	fprintf(stderr,"normalize_probs\n");
	print_probs(stderr, P);
	fprintf(stderr,"normalize_probs done\n");
#endif
}


/* Make a primitive PSSM from a sequence with up to M mismatches.
   Score "match" for match and "mismatch" for mismatch
   Note that alphsize is the length NOT including the wildcard (last
   letter in alphabet, numbered alphsize).
*/
PSSM string_to_pssm(ubyte_t *seq, int len, int alphsize, float match, float mismatch,
			float wcscore) {
  PSSM mat;
  int i, j, nScores;
  int *scores, *base;

  nScores = len*(alphsize+1);

  /* initialize scores to mismatch */
  base = scores = (int *)malloc(nScores*sizeof(int));
  for (i=0; i<len; ++i) for (j=0; j<=alphsize; ++j) *(scores++) = (int)(1000 * mismatch);

  /* Score for wildcards in genome */
  scores = base;
  for (i=0; i<len; ++i) {
    scores[alphsize] = wcscore;
    scores += alphsize+1;
  }

  /* Make matches score 0 */
  scores = base;
  for (i=0; i<len; ++i) {
    /* Wildcards in query*/
    if (seq[i]==alphsize) for (j=0; j<=alphsize; ++j) scores[j] = wcscore;
    else scores[seq[i]] = match;
    scores += alphsize+1;
  }

  /* This function makes a matrix */
  mat = init_matrix_score(0, len, alphsize+1, base, nScores, -0.5);

  return mat;
}


/* Allocates a Markov chain in a Probs structure
   len = (alphsize+1)^order
   It means that the total array length is len*(alphsize+1)
*/
Probs *alloc_markov_chain(int order, int alphsize) {
  int i;
  int *powers;
  Probs *P;

  /* powers array holds powers of alphsize+1 from 0 to order+1 */
  powers = (int *)malloc((order+2)*sizeof(int));
  powers[0]=1;
  for (i=1; i<order+2; ++i) powers[i]=(alphsize+1)*powers[i-1];

  P = alloc_probs(powers[order], alphsize);
  P->powers = powers;
  P->order = order;
  return P;
}


/*
 * The translation from ascii quality score to probability is stored in a
 * look-up table
 *
 * Solexa style qual scores (now they use Phred -?)
 */
float *solexa_ascii_quality_scores() {
	float *errorprob = calloc(128,sizeof(float));
	int low = 33;
	int high = 126;
	int i;
	for (i=0; i<low; ++i) errorprob[i] = 0.75;
	for (i=low; i<=high; ++i) {
		errorprob[i] = 1./(1.+pow(10.,0.1*(i-64)));
		if (errorprob[i]>0.75) errorprob[i] = 0.75;
	}
	for (i=high+1; i<128; ++i) errorprob[i] = 0.;
	return errorprob;
}


/*
 * Phred style qual scores
 * base is normally 33 or 64
 */
float *phred_ascii_quality_scores(int base) {
	float *errorprob = calloc(128,sizeof(float));
	int low = base;
	int high = 126;
	int i;
	for (i=0; i<low; ++i) errorprob[i] = 0.75;
	for (i=low; i<=high; ++i) {
		errorprob[i] = pow(10.,-0.1*(i-base));
		if (errorprob[i]>0.75) errorprob[i] = 0.75;
	}
	for (i=high+1; i<128; ++i) errorprob[i] = 0.;
	return errorprob;
}

/*
 * Load the error model from a file with the following format:
 *
 * 
 * # comments
 * # Each line has qual score and base followed by 4 scores
 * # for A, C, G and T. Fields separated by blanks
 * # Like this
 * 0 A 0. 0. 0. 0.
 * 0 C 0. 0. 0. 0.
 * :
 * 3 A 1.3 -2.1 -3.1 -2.0
 * 3 C -1.9 2.0 -2.1 -1.8
 * :
 * :
 * 50 T ..
 */
void load_error_model(float *table, const char *filename) {
    FILE *fp = fopen(filename, "r");
    char line[1024];

    int index = -1;
    char base;
    int base_index = -1;

    if (!fp) {
        fprintf(stderr, "Error model file %s not found.\n", filename);
        exit(1);
    }
    
    while(fgets(line, 1024, fp) != NULL)
    {
        char *pch;
        int counter = 0;

        if (line[0] == '#')
            continue;
        
        pch = strtok(line, " \t");
        while (pch != NULL)
        {
            //fprintf(stderr, "counter: %d pch: %s\n", counter, pch);
            switch(counter) {
                case 0:
                    index = atoi(pch);
                    break;
                case 1:
                    base = pch[0];
                    base_index = nst_nt4_table[(int)base];
                    break;
                default:
                    //fprintf(stderr, "index: %d base_index: %d counter-2: %d total_index: %d\n", index, base_index, counter-2, total_index);
                    table[index * 16 + base_index * 4 + counter - 2] = atof(pch);
            }

            pch = strtok(NULL, " \t");
            counter++;
        }
    }

    fclose(fp);
}

/*
 * Read qual values in a file with two columns: ASCII code and prob
 * Anything else will fail. Sorry. E.g.:
 * B 0.23456
 * C 0.20743
 * D 0.18394
 * ...
 * Print out the result in compressed format
 */
float *read_ascii_quality_scores(char *filename) {
	float *errorprob, p, last;
    double pd;
	int i, n;
	char c;
	FILE *fp = fopen(filename,"r");
	if (!fp) {
		fprintf(stderr,"ERROR: Couldn't open file %s for reading qualities\n",filename);
		exit(101);
	}
	errorprob = calloc(128,sizeof(float));
	for (n=0; n<128; ++n) errorprob[n] = -1.;
	n=0;
	fprintf(stderr,"# qual: ");
	while ( fscanf(fp,"%1s %lf",&c,&pd) != EOF ) {
        p = (float)pd;
		i = (int)c;
		if (i<0 || i>=128 || p>1. || p<0.) {
			fprintf(stderr,"ERROR: These quals are strange %c %f\n",i,p);
			exit(102);
		}
		errorprob[i] = p;
		if (n%15==15) fprintf(stderr,"\n# qual: ");
		else fprintf(stderr," ");
		fprintf(stderr,"%c%5.4f",i,errorprob[i]);
		++n;
	}
	fprintf(stderr,"\n");
	last = 0.75;
	for (n=0; n<128; ++n) {
		if (errorprob[n] < 0.) errorprob[n]=last;
		else last = errorprob[n];
	}
	return errorprob;
}



/* 
 * Make 0'th order PSSM given the scores in the error model.
 *
 * Each base quality will have a set of scores associated
 * with it, and each of these will be added to the matrix.
 */
PSSM error_model_to_pssm(PSSM mat, ubyte_t *seq, ubyte_t *qual, int len, int alphsize,
        const float *error_model, int fastq_base) {

    int i, k, q;

    for (i = 0; i < len; i++) {
        /* The error prob from the qual ascii code */
        q = (int)qual[i];
        if (q<0 || q>=128) {
            fprintf(stderr,"Weird qual: %d\n%s\n%s",q,seq,qual);
        }
        /*fprintf(stderr,"q:  %d      i:  %d \n",q-33,i);*/
        for (k = 0; k < alphsize; k++) {
            mat->scores[mat->offsets[i] + k] = (int) (1000 * error_model[16 *
                    (q - 33) + 4 * seq[i] + k]);
        }
    }

    return mat;
}

/* Turn quality scores into probabilities for each letter
   and return in an array of size len*(alphsize+1).

   alphsize is the actual alphabet size (4 for DNA).
   A wildcard char of alphasize is assumed
*/
Probs *qual_to_probs(ubyte_t *seq, ubyte_t *qual, int len, int alphsize,
		float *qualprobs) {
	int i, j, q;
	Probs *P;
	float *prob, p1, p2;

	P = alloc_probs(len, alphsize);
	prob = P->p;
	for (i = 0; i < len; ++i) {
		/* The error prob from the qual ascii code */
		q = (int)qual[i];
		if (q<0 || q>=128) {
			fprintf(stderr,"Weird qual: %d\n%s\n%s",q,seq,qual);
		}
		p1 = qualprobs[(int)qual[i]];
		/* The prob of the letter called */
		p2 = p1 / 3.;
		for (j = 0; j < seq[i] && j < alphsize; ++j) prob[j] = p2;
		prob[j] = 1. - p1;
		for (j = j + 1; j < alphsize; ++j) prob[j] = p2;
		/* If query is a wildcard, the probs are uniform */
		if (seq[i] == alphsize) for (j = 0; j < alphsize; ++j)
			prob[j] = 1. / alphsize;
		/* Wildcards in genome set to prob 1/alphsize */
		prob[alphsize] = 1. / alphsize;
		prob += alphsize + 1;
	}

#ifdef DEBUG
	fprintf(stderr,"qual_to_probs\n");
	print_probs(stderr, P);
	fprintf(stderr,"qual_to_probs done\n");
#endif

	return P;
}



/* Modify an array of probabilities with the prior prob of an
   error (in the reference genome) or a SNP (psnp) according to
   this formula (for each position):
     P(base b) = (1-psnp/[1-q(b)])P_0(b) + psnp q(b) sum_c P_0(c)/[1-q(c)]
   P_0(c) is the probability before SNPs
   q is the base frequency - it is assumed that the base is
   randomly drawn from this dist if there is a SNP or error.
   define f(c) = psnp /[1-q(c)], then
     P(b) = ( 1 - f(b) ) P_0(b) + q(b) sum_c P_0(c) f(c)
*/
void snp_probs(Probs *P, float *q, float psnp) {
  int i, b;
  float bg[256], f[256], p, sum, *prob;

  prob = P->p;

  if (!q) {
    p = 1.0/P->alphsize;
    for (b=0; b<P->alphsize; ++b) bg[b]=p;
    q=bg;
  }
  for (b=0; b<P->alphsize; ++b) f[b]=psnp/(1.-q[b]);

  for (i=0; i<P->len; ++i) {
    for (sum = 0., b=0; b<P->alphsize; ++b) sum += prob[b]*f[b];
    for (b=0; b<P->alphsize; ++b) prob[b] = (1.-f[b])*prob[b]+q[b]*sum;
    prob += P->alphsize+1;
  }
#ifdef DEBUG
  fprintf(stderr,"snp_probs\n");
  print_probs(stderr, P);
  fprintf(stderr,"snp_probs done\n");
#endif
}



/* Make an n'th order PSSM from positional base probabilities
   background is the n-th order Markov Chain in logarithmic form:
         log(x_i|x_{i-n},x_{i-n+1},...,x_{i-1})
   It is assumed that the background is set to log(1/alphsize)
   for any word with wildcards.

   If order==0 && background==NULL a uniform background is assumed

   A matrix of len+order is initialized.
   All scores in the first order positions are set to 0.

   alphsize is the actual alphabet size (4 for DNA).
   A wildcard char of alphasize is assumed
*/
PSSM prob_to_pssm(Probs *P, Probs *mc) {
  PSSM mat;
  int i, j, N, order=0;
  int *scores;
  float *prob, pflat;
  float *background, bg[256];

  if (mc) order = mc->order;

  mat = init_matrix(order, P->len+order, P->alphsize+1);
  scores = mat->scores;

  /* Initialize beginning of matrix */
  if (order) {
    if (!mc) {
      fprintf(stderr,"No MC chain supplied for n-th order matrix\n");
      return NULL;
    }
    for (i=0; i<mat->offsets[order]; ++i) scores[i] = 0;
  }
  if (!mc) {
    pflat = log2f(1.0/P->alphsize);
    for (i=0; i<P->alphsize; ++i) bg[i]=pflat;
    background=bg;
  }
  else background = mc->p;

  /* Filling in scores for each letter */
  prob=P->p;
  for (i=0; i<P->len; ++i) {
    /* Pointer to the beginning of the score array */
    scores = mat->scores+mat->offsets[i+order];
    N = mat->offsets[i+order+1] - mat->offsets[i+order];
    for (j=0; j<N; ++j) {
      scores[j]= (int)(1000 * (log2f(prob[j%(P->alphsize+1)]) - background[j]));
    }
    /* Pointer to the appropriate probabilities */
    prob += P->alphsize+1;
  }

  return mat;
}


/* Calculate the maximum possible score for a PSSM */
int highest_scores(PSSM mat) {
  int i, j, N, max;
  int order = mat->order;
  int hscore;
  int *scores = mat->scores;

  hscore=0.;
  for (i=order; i<mat->length; ++i) {
    /* Pointer to the beginning of the score array */
    scores = mat->scores+mat->offsets[i];
    /* Find max score */
    N = mat->offsets[i+1] - mat->offsets[i];
    max=0;
    for (j=1; j<N; ++j) if (scores[j]>=scores[max]) max=j;
    hscore += scores[max];
  }
  return hscore;
}




/* Set matrix threshold such that at most M mismatches are allowed.
   Finds the M cheapest non-consensus scores and calculates threshold.
*/
int mismatch_threshold(PSSM mat, int M) {
	int i, j, k, N, cheapest[MAXPSSMSIZE], max[2];
	int order = mat->order;
	int scorediff[MAXPSSMSIZE], t, hscore;
	int *scores = mat->scores;
	const int infty = 1.e-100;

	hscore = 0.;
	for (i = 0; i < mat->length; ++i)
		scorediff[i] = infty;
	for (i = order; i < mat->length; ++i) {
		/* Pointer to the beginning of the score array */
		scores = mat->scores + mat->offsets[i];
		/* Find max score and second highest */
		N = mat->offsets[i + 1] - mat->offsets[i];
		if (scores[0] > scores[1]) {
			max[0] = 0;
			max[1] = 1;
		}
		else {
			max[0] = 1;
			max[1] = 0;
		}
		for (j = 2; j < N; ++j) {
			if (scores[j] >= scores[max[1]]) {
				if (scores[j] >= scores[max[0]]) {
					max[1] = max[0];
					max[0] = j;
				}
				else max[1] = j;
			}
		}
		hscore += scores[max[0]];
		scorediff[i] = scores[max[0]] - scores[max[1]];
		k = i;
		cheapest[i] = i;
		while (k > order && scorediff[cheapest[k]] < scorediff[cheapest[k - 1]]) {
			j = cheapest[k];
			cheapest[k] = cheapest[k - 1];
			cheapest[k - 1] = j;
			--k;
		}
	}

	t = hscore;
	for (k = 0; k < M; ++k)
		t -= scorediff[cheapest[k + order]];
	t -= 0.5 * scorediff[cheapest[k + order]];

#ifdef DEBUG
	for (i=order; i<mat->length; ++i) fprintf(stderr,"%f %d %f\n",
			scorediff[i],cheapest[i],scorediff[cheapest[i]]);
	fprintf(stderr,"mismatch_threshold, max score %f, threshold: %f\n",hscore,t);
#endif

	return t;
}


/* Makes a Markov chain from an array of counts
   Log probabilities are stored in one IndexType array following
   the (order+1) word numbers (words of letters from 0 to
   AlphLen).
*/
static void Counts2markov_chain(Probs *mc) {
  int i, j, k, wild, al, lett[256], *counts;
  float *p, sum, pflat;

  pflat = (float)log2f((float)1.0f/mc->alphsize);

  /* Initialize lett array, which keeps track of the word corresponding
     to any given number (spedometer)
  */
  for (k=0;k<=mc->order;++k) lett[k]=0;

  al = mc->alphsize;
  p = mc->p;
  counts = mc->counts;
  for (i=0; i< mc->len; ++i) {
    /* Check for wildcards anywhere in word */
    for (k=0, wild=0; k<mc->order; ++k) if (lett[k]==al) { wild=1; break; }

    /* All probs equal to pflat for wildcards (perhaps not best choice) */
    if (wild) for (j=0; j<al; ++j) *(p++) = pflat;
    else {
      sum = 0.;
      for (j=0; j<al; ++j) sum += counts[j];
      for (j=0; j<al; ++j) *(p++) = log2f(counts[j]/sum) ;
    }
    /* Prob of a wildcard in the genome is always pflat */
    *(p++) = pflat;
    counts += al+1;

    /* Keep track of the word */
    for (k=mc->order-1; k>=0; --k) {
      lett[k]+=1;
      if (lett[k]>al) lett[k]=0;
      else break;
    }
  }
#ifdef DEBUG
  fprintf(stderr,"Counts2markov_chain\n");
  print_probs(stderr, mc);
  fprintf(stderr,"Counts2markov_chain done\n");
#endif
}


Probs *markov_chain(bwtint_t *counts, int alphlen) {
  Probs *mc = alloc_markov_chain(0, alphlen);
  int i;
  mc->counts = malloc(alphlen * sizeof(int));
  for (i = 0; i < alphlen; i++) { 
    mc->counts[i] = counts[i+1] - counts[i];
    //mc->counts[i] = counts[alphlen] / alphlen;
  }
  Counts2markov_chain(mc);
  return mc;
}

void set_thresholds(PSSM mat, const gap_opt_t *opt)
 {
    int i;
    ubyte_t j, k;
    int biggest_drop = -INT_MAX;
    int best_score, total_best_score = 0.0;
    int seed_best_score = -1.0;
    int drop;

    for (i = 0; i < get_length_fast(mat); i++) {
        best_score = -INT_MAX;

        for (j = 0; j < 4; j++) {
            int score1 = get_score_fast(mat, &j, i);

            if (score1 > best_score)
                best_score = score1;

            for (k = j + 1; k < 4; k++) {
                int score2 = get_score_fast(mat, &k, i);
                drop = score1 - score2;
                if (drop < 0)
                    drop = -drop;

                if (drop > biggest_drop)
                    biggest_drop = drop;
            }
        }

        if (i == opt->seed_len)
            seed_best_score = total_best_score;

        total_best_score += best_score;
    }

    biggest_drop = 15000;

    if (opt->pssm_ratio > 0.0)
        calc_and_set_reverse_thresholds(mat, 1, get_length(mat), total_best_score - opt->pssm_ratio * biggest_drop);
    else
        calc_and_set_reverse_thresholds(mat, 1, get_length(mat), opt->threshold);

    if (seed_best_score > 0 && opt->pssm_seed_ratio > 0) {
        //set_length(mat, opt->seed_len);
        calc_and_set_reverse_thresholds(mat, get_length(mat) - opt->seed_len, get_length(mat), seed_best_score - opt->pssm_seed_ratio * biggest_drop);
        //set_length(mat, len);
    }
 }


/* Make a matrix from a sequence
 *
 * scoretype determines the type:
 *
 * 0: Do NOT use quality scores. Use just sequence with match, mismatch
 *    and wildcard scores
 * MATRIX_TYPE: A matrix has been read and just needs normalizing etc
 * 	  (zero order is assumed)
 * 	  psnp is used as pseudocount!
 * Otherwise: Use qual. scores. Uses a uniform letter prob for SNPs
 *    (NULL in snp_probs)
 *    The mc chain background may be NULL (means 0th order uniform)
 *
*/
int sequence_to_pssm(bwa_seq_t *s, int alphsize, float psnp, Probs *mc, float sc_match,
		float sc_mismatch, float sc_wild, int scoretype, float *qualprobs,const gap_opt_t *opt)
 {
	int nf=0, nr=0;
    int i;
	Probs *P;

	if (scoretype==0) {
		s->mat = string_to_pssm(s->seq, s->len, alphsize, sc_match, sc_mismatch, sc_wild);
		//if (s->rseq) s->revmat = string_to_pssm(s->rseq, s->len, alphsize, sc_match, sc_mismatch, sc_wild);
	}
	else if (scoretype==MATRIX_TYPE) {
		P = alloc_probs(s->len,alphsize);
		free(P->p);
		P->p = (float *)(s->qual);
		s->qual = NULL;
		normalize_probs(P,psnp);
		s->mat = prob_to_pssm(P, mc);
		free_probs(P);
	}
	else {
		// Remove Ns in beginning of sequence:
        // while (nf < s->len && s->seq[nf]>=alphsize) ++nf;
        // Causes segfaults, and it doesn't really hurt anything.
        // Besides, with the BW implementation, the beginning of the
        // sequence is really the end

        P = qual_to_probs(s->seq+nf, s->rqual+nf, s->len-nf, alphsize, qualprobs);
        snp_probs(P, NULL, psnp);

        if (opt->parclip) {
            for (i = 0; i < P->len; i++) {
                // set the probability of C and T to be equal to their average
                // this isn't strictly correct, but oh well
                // float avg = P->p[i*P->alphsize + 1] + P->p[i*P->alphsize + 3];
                float big_val = P->p[i*P->alphsize + 1] > P->p[i*P->alphsize + 3] ? P->p[i*P->alphsize + 1] : P->p[i*P->alphsize + 3];
                float small_val = P->p[i*P->alphsize + 1] > P->p[i*P->alphsize + 3] ? P->p[i*P->alphsize + 1] : P->p[i*P->alphsize + 3];
                // fprintf(stderr, "big_val: %f small_val: %f\n", big_val, small_val);

                P->p[i*P->alphsize + 3] = small_val;
                P->p[i*P->alphsize + 1] = big_val;
            }
        }

        s->mat = prob_to_pssm(P, mc);
        free_probs(P);

        if (opt->use_error_model) 
            error_model_to_pssm(s->mat, s->seq+nf, s->rqual+nf, s->len-nf, alphsize, opt->error_lookup, opt->fastq_base);

        /*
        if (debug)
            fprintf(stderr, "qual: %s\n", s->rqual);

        for (i = 0; i < s->len; i++) {
            if (debug)
                fprintf(stderr, "%d qual: %d seq: %d || ", i, s->rqual[i] - 33, s->seq[i]);

            for (j = 0; j < 4; j++) {
                if (debug)
                    fprintf(stderr, "%d ", get_score_fast(s->mat, &j, i) / 100);
            }

            if (debug)
                fprintf(stderr, "\n");
        }
        */
	}

    //set_thresholds(s->mat, opt);
	return nf+nr;
}

