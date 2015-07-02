#ifndef TILEHMM_LIB_H_
#define TILEHMM_LIB_H_



#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP _log_sum(SEXP x, SEXP y, SEXP b);
SEXP _forward(SEXP hmm, SEXP obs, SEXP nstates, SEXP obs_len);
SEXP _backward(SEXP hmm, SEXP obs, SEXP nstates, SEXP obs_len);
SEXP _baumWelch_trans(SEXP hmm, SEXP obs, SEXP alpha, SEXP beta, SEXP trans_prior, SEXP init_prior);
SEXP _viterbi(SEXP hmm, SEXP obs);
SEXP __calcTau(SEXP gamma, SEXP obs, SEXP hmm);
SEXP __get_xi(SEXP hmm, SEXP alpha, SEXP beta, SEXP trans, SEXP obs);
SEXP __get_gamma(SEXP xi_t, SEXP alpha);
SEXP __get_obs_prob(SEXP hmm, double obs);
SEXP _get_obs_prob(SEXP hmm, SEXP obs);
#endif /*TILEHMM_LIB_H_*/
