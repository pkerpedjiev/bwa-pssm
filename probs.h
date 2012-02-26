#ifndef PROBS_H
#define PROBS_H 1

/* This struct holds probabilities or a Markov Model */
typedef struct {
  int len;
  int alphsize;
  int order;
  int *powers;
  int *counts;
  float *p;
} Probs;


#endif
