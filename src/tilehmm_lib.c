#include "tilehmm_lib.h"

/* maximum of a double array */
double* max(double *array, int length){
	double *max_entry;
	max_entry = malloc(2 * sizeof(double));
	max_entry[0] = array[0];
	max_entry[1] = 0;
	for(int i=1; i < length; i++){
		if(array[i] > max_entry[0]){
			max_entry[0] = array[i];
			max_entry[1] = (double) i;
		}			
	} 
		
	return max_entry;
} 


/* calculate log(x + y) from log(x) and log(y) */
SEXP _log_sum(SEXP x, SEXP y, SEXP b){
		double *c_x, *c_y, log_term, min_arg, max_arg;
		int *c_b;
		SEXP ret;
		
		PROTECT(ret = allocVector(REALSXP,1));
		c_x = REAL(x);
		c_y = REAL(y);
		c_b = INTEGER(b);
		
		if(*c_x == -INFINITY){
			REAL(ret)[0] = *c_y;
		}
		else if(*c_y == -INFINITY){
			REAL(ret)[0] = *c_x;
		}
		else{
			if(*c_x <= *c_y){ 
				min_arg = *c_x;
				max_arg = *c_y;
			}
			else{
				min_arg = *c_y;
				max_arg = *c_x;
			} 
			
			if(*c_b == 10) log_term = log10(1 + pow(10,min_arg - max_arg));
			else if(*c_b == 0) log_term = log(1 + exp(min_arg - max_arg));
			else log_term = log(1 + pow(*c_b,min_arg - max_arg))/log(*c_b);
			
			REAL(ret)[0] = max_arg + log_term;
		}
		UNPROTECT(1);
		return(ret);
}

/* calculate forward variables 
 * only for HMMs derived from contHMM */
SEXP _forward(SEXP hmm, SEXP obs, SEXP nstates, SEXP obs_len){
	int c_rows = INTEGER(nstates)[0], c_cols = INTEGER(obs_len)[0];
	double *c_alpha, *c_init, *c_trans;
	SEXP alpha, init_slot, log_trans, trans_slot, obs_prob, tmp_alpha, tmp_alpha2;
	SEXP base;
	PROTECT(base = allocVector(INTSXP,1));
	INTEGER(base)[0] = 0;
	
	/* get initial probability distribution for states */
	PROTECT(init_slot = allocVector(STRSXP,1));
	SET_STRING_ELT(init_slot,0,mkChar("init"));
	c_init = REAL(R_do_slot(hmm,init_slot));
	UNPROTECT(1); // init_slot
	
	/* get log transformed transition probabilities */
	PROTECT(trans_slot = allocVector(STRSXP,1));
	SET_STRING_ELT(trans_slot,0,mkChar("transition.matrix"));	
	c_trans = REAL(R_do_slot(hmm,trans_slot));
	UNPROTECT(1); // trans_slot
	PROTECT(log_trans = allocMatrix(REALSXP,c_rows,c_rows));
	for(int i=0; i < c_rows * c_rows; i++){
		if(c_trans[i] == 0){
			REAL(log_trans)[i] = -INFINITY;
		}
		else{
			REAL(log_trans)[i] =  log(c_trans[i]);
		}
	}
	
	/* matrix for forward variables */
	PROTECT(alpha = allocMatrix(REALSXP,c_rows,c_cols));
	c_alpha = REAL(alpha);
	
	/* initialise */
	PROTECT(obs_prob = __get_obs_prob(hmm,REAL(obs)[0]));
	for(int i=0; i < c_rows; i++){
		if(c_init[i] == 0){
			c_alpha[i] = -INFINITY;
		}
		else{
			c_alpha[i] = log(c_init[i]) + REAL(obs_prob)[i];
		}
	}
	UNPROTECT(1);  // obs_prob
	
	/* induction */
	PROTECT(tmp_alpha = allocVector(REALSXP,1));
	PROTECT(tmp_alpha2 = allocVector(REALSXP,1));
	for(int t = 1; t < c_cols; t++){
		PROTECT(obs_prob = __get_obs_prob(hmm,REAL(obs)[t]));
		for(int j = 0; j < c_rows; j++){
			/* i = 0 */
			REAL(tmp_alpha)[0] = c_alpha[(t-1)*c_rows] + REAL(log_trans)[j*c_rows];
			for(int i = 1; i < c_rows; i++){
				REAL(tmp_alpha2)[0] = c_alpha[i + (t-1)*c_rows] + REAL(log_trans)[i + j*c_rows];
				REAL(tmp_alpha)[0] = REAL(_log_sum(tmp_alpha, tmp_alpha2, base))[0];
			}
			c_alpha[j + t*c_rows] = REAL(tmp_alpha)[0] + REAL(obs_prob)[j];
		}
		UNPROTECT(1); //obs_prob
	}
	
	UNPROTECT(5); // log_trans, alpha, tmp_alpha, tmp_alpha2, base
	return alpha;
}

/* calculate backward variables 
 * only for HMMs derived from contHMM */
SEXP _backward(SEXP hmm, SEXP obs, SEXP nstates, SEXP obs_len){
	int c_rows = INTEGER(nstates)[0], c_cols = INTEGER(obs_len)[0];
	double *c_beta, *c_trans;
	SEXP beta, log_trans, trans_slot, obs_prob, tmp_beta, tmp_beta2;
	SEXP base;
	PROTECT(base = allocVector(INTSXP,1));
	INTEGER(base)[0] = 0;
	
	/* get log transformed transition probabilities */
	PROTECT(trans_slot = allocVector(STRSXP,1));
	SET_STRING_ELT(trans_slot,0,mkChar("transition.matrix"));	
	c_trans = REAL(R_do_slot(hmm,trans_slot));
	UNPROTECT(1); // trans_slot
	PROTECT(log_trans = allocMatrix(REALSXP,c_rows,c_rows));
	for(int i=0; i < c_rows * c_rows; i++){
		if(c_trans[i] == 0){
			REAL(log_trans)[i] = -INFINITY;
		}
		else{
			REAL(log_trans)[i] =  log(c_trans[i]);
		}
	}
	
	/* matrix for backward variables */
	PROTECT(beta = allocMatrix(REALSXP,c_rows,c_cols));
	c_beta = REAL(beta);
	
	/* initialise */
	for(int i=0; i < c_rows; i++){
		c_beta[i + (c_cols-1)*c_rows] = 0;
	}
	
	/* induction */
	PROTECT(tmp_beta = allocVector(REALSXP,1));
	PROTECT(tmp_beta2 = allocVector(REALSXP,1));
	for(int t = c_cols-2; t >= 0; t--){
		PROTECT(obs_prob = __get_obs_prob(hmm,REAL(obs)[t+1]));
		for(int i = 0; i < c_rows; i++){
			REAL(tmp_beta)[0] = c_beta[(t+1)*c_rows] + REAL(obs_prob)[0] + REAL(log_trans)[i];
			for(int j = 1; j < c_rows; j++){
				REAL(tmp_beta2)[0] = c_beta[j + (t+1)*c_rows] + REAL(obs_prob)[j] + REAL(log_trans)[i + j*c_rows];
				REAL(tmp_beta)[0] = REAL(_log_sum(tmp_beta,tmp_beta2,base))[0];
			}
			c_beta[i + t*c_rows] = REAL(tmp_beta)[0];
		}
		UNPROTECT(1); //obs_prob
	}
	
	UNPROTECT(5); // log_trans, beta, tmp_beta, tmp_beta2, base
	return beta;
}

/* Computes xi variables for one observation sequence. 
 * Returns a list of matrices. 
 * Note that this only calculates xi for t = 1, ..., length(obs)-1 */ 
SEXP __get_xi(SEXP hmm, SEXP alpha, SEXP beta, SEXP trans, SEXP obs){
	SEXP xi_t, dim_attr, xi_ij, xi_denom, xi_num, obs_prob;
	SEXP base, alpha_tmp;
	int c_dim;
	
	PROTECT(base = allocVector(INTSXP,1));
	INTEGER(base)[0] = 0;
	
	/* get matrix dimensions */
	PROTECT(dim_attr = allocVector(STRSXP,1));
	SET_STRING_ELT(dim_attr,0,mkChar("dim"));
	c_dim = INTEGER(getAttrib(trans,dim_attr))[0];
	UNPROTECT(1); /* dim_attr */
	
	/* initialise list */
	PROTECT(xi_t = allocVector(VECSXP, length(obs)-1));
	PROTECT(xi_denom = allocVector(REALSXP,1));
	PROTECT(xi_num = allocVector(REALSXP,1));
	PROTECT(alpha_tmp = allocVector(REALSXP,1));
	
	/* compute xi */
	for(int t=0; t < length(obs)-1; t++){
		PROTECT(xi_ij = allocMatrix(REALSXP,c_dim,c_dim));
		PROTECT(obs_prob = __get_obs_prob(hmm,REAL(obs)[t+1]));		
		for(int i=0; i < c_dim; i++){
			for(int j=0; j < c_dim; j++){
				REAL(xi_num)[0] = REAL(alpha)[i + t*c_dim] + REAL(trans)[i + j*c_dim] + REAL(obs_prob)[j] + REAL(beta)[j + (t+1)*c_dim];
				REAL(xi_ij)[i + j*c_dim] = REAL(xi_num)[0]; 
			}
		}
		
		/* normalise xi_ij*/
		REAL(xi_denom)[0] = REAL(alpha)[(length(obs)-1)*c_dim];
		for(int i = 1; i < c_dim; i++){
			REAL(alpha_tmp)[0] = REAL(alpha)[i + (length(obs)-1)*c_dim];
			REAL(xi_denom)[0] = REAL(_log_sum(xi_denom,alpha_tmp,base))[0];
		}
		for(int i=0; i < c_dim; i++){
			for(int j=0; j < c_dim; j++){
				REAL(xi_ij)[i + j*c_dim] = REAL(xi_ij)[i + j*c_dim] - REAL(xi_denom)[0];
			}
		}
		SET_VECTOR_ELT(xi_t,t,xi_ij);
		UNPROTECT(2); /* xi_ij, obs_prob */
	}
	
	UNPROTECT(5); /* xi_t, xi_denom, xi_num, alpha_tmp, base */
	return xi_t;
}

/* compute gamma from precomputed xi and alpha (for t = T) */
SEXP __get_gamma(SEXP xi_t, SEXP alpha){
	SEXP gamma_t, xi_ij, gamma_denom, dim_attr, dim, base, tmp, tmp2;
	int c_len, c_dim;
	
	PROTECT(base = allocVector(INTSXP,1));
	INTEGER(base)[0] = 0;
	PROTECT(tmp = allocVector(REALSXP,1));
	PROTECT(tmp2 = allocVector(REALSXP,1));
	
	/* get matrix dimension */
	PROTECT(dim_attr = allocVector(STRSXP,1));
	SET_STRING_ELT(dim_attr,0,mkChar("dim"));
	PROTECT(dim = getAttrib(alpha,dim_attr));
	c_dim = INTEGER(dim)[0];
	c_len = INTEGER(dim)[1];
	UNPROTECT(2); /* dim_attr, dim */
	
	/* compute gamma for t= 1, ..., T-1 */ 
	PROTECT(gamma_t = allocMatrix(REALSXP,c_dim,c_len));
	for(int t=0; t < c_len-1; t++){
		xi_ij = VECTOR_ELT(xi_t,t);
		for(int i=0; i < c_dim; i++){
			REAL(gamma_t)[i + t*c_dim] = -INFINITY;
			for(int j=0; j < c_dim; j++){
				REAL(tmp)[0] = REAL(gamma_t)[i + t*c_dim];
				REAL(tmp2)[0] = REAL(xi_ij)[i + j*c_dim];
				REAL(gamma_t)[i + t*c_dim] = REAL(_log_sum(tmp,tmp2,base))[0];
			}
		}
	}
	/* compute gamma for t = T */
	PROTECT(gamma_denom = allocVector(REALSXP,1));
	REAL(gamma_denom)[0] = REAL(alpha)[(c_len-1)*c_dim];
	for(int i=1; i < c_dim; i++){
//		REAL(tmp)[0] = REAL(gamma_denom)[0];
		REAL(tmp2)[0] = REAL(alpha)[i + (c_len-1)*c_dim];
		REAL(gamma_denom)[0] =  REAL(_log_sum(gamma_denom, tmp2, base))[0];
	}
	for(int i=0; i < c_dim; i++){
		REAL(gamma_t)[i + (c_len-1)*c_dim] = REAL(alpha)[i + (c_len-1)*c_dim] - REAL(gamma_denom)[0];
	}	
	
	UNPROTECT(5); /* gamma_t, base, gamma_denom, tmp, tmp2 */
	return gamma_t;
}

/* Calculate new estimates for transition and initial probabilities 
 * obs, alpha, and beta are lists with one element for each observation sequence.
 * prior is a matrix with priors for transition probability distributions
 * Returns list with components 'pi', 'transition', and 'gamma' */
SEXP _baumWelch_trans(SEXP hmm, SEXP obs, SEXP alpha, SEXP beta, SEXP trans_prior, SEXP init_prior){
	SEXP ret, pi, trans, log_trans, gamma, xi, gamma_t, xi_t, gamma_sum_t, xi_sum_t;
	SEXP base, dim_attr, trans_slot, tmp1, tmp2, xi_d, gamma_d, ret_names, prior_sum;
	int c_D, c_T, c_dim;
	double c_sum;

	PROTECT(base = allocVector(INTSXP,1));
	INTEGER(base)[0] = 0;
	PROTECT(tmp1 = allocVector(REALSXP,1));
	PROTECT(tmp2 = allocVector(REALSXP,1));
	
	/* get current transition probabilities */
	PROTECT(trans_slot = allocVector(STRSXP,1));
	SET_STRING_ELT(trans_slot,0,mkChar("transition.matrix"));	
	trans = R_do_slot(hmm,trans_slot);
	UNPROTECT(1); /* trans_slot */
	
	/* get dimensions */
	c_D = length(obs);
	PROTECT(dim_attr = allocVector(STRSXP,1));
	SET_STRING_ELT(dim_attr,0,mkChar("dim"));
	c_dim = INTEGER(getAttrib(trans,dim_attr))[0];
	UNPROTECT(1); /* dim_attr */
	
	/* log transform transition probabilities*/
	PROTECT(log_trans = allocMatrix(REALSXP,c_dim,c_dim));
	for(int i=0; i < c_dim * c_dim; i++){
		REAL(log_trans)[i] =  REAL(trans)[i];
		if(REAL(log_trans)[i] == 0){
			REAL(log_trans)[i] = -INFINITY;
		}
		else{
			REAL(log_trans)[i] = log(REAL(log_trans)[i]);
		}
	}	
	
	/* calculate sum of prior distributions for each state
	 * note that they may be scaled differently and are not guaranteed to sum to 1 */
	PROTECT(prior_sum = allocVector(REALSXP, c_dim));
	for(int i=0; i < c_dim; i++){
		c_sum = 0;
		for(int j=0; j < c_dim; j++){
			c_sum += REAL(trans_prior)[i + j*c_dim]; 
		}
		REAL(prior_sum)[i] = c_sum;
	}
	/* calculate sum of prior initial state probabilities */
	c_sum = 0;
	for(int i=0; i < c_dim; i++){
		c_sum += REAL(init_prior)[i];
	}
	
	/* create and populate lists for gamma and xi */
	PROTECT(xi = allocVector(VECSXP, c_D));
	PROTECT(gamma = allocVector(VECSXP, c_D));
	
	for(int d=0; d < c_D; d++){
		SET_VECTOR_ELT(xi,d,__get_xi(hmm,VECTOR_ELT(alpha,d),VECTOR_ELT(beta,d), log_trans, VECTOR_ELT(obs,d)));
		SET_VECTOR_ELT(gamma, d, __get_gamma(VECTOR_ELT(xi, d),VECTOR_ELT(alpha,d)));
	}

	/* sum xi and gamma over t for all observation sequences */
	PROTECT(xi_sum_t = allocVector(VECSXP,c_D));
	PROTECT(gamma_sum_t = allocVector(VECSXP,c_D));
	PROTECT(pi = allocVector(REALSXP,c_dim));
	
	for(int d=0; d < c_D; d++){
		PROTECT(xi_t = allocMatrix(REALSXP,c_dim,c_dim));
		PROTECT(gamma_t = allocVector(REALSXP,c_dim));
		c_T = length(VECTOR_ELT(obs,d));
			
		/* initialise new elements*/
		for(int i=0; i < c_dim; i++){
			REAL(gamma_t)[i] = REAL(VECTOR_ELT(gamma,d))[i]; 
			for(int j=0; j < c_dim; j++){
				REAL(xi_t)[i + j*c_dim] = REAL(VECTOR_ELT(VECTOR_ELT(xi,d),0))[i + j*c_dim];
			}
		}
		/* compute sum */
		for(int t=1; t < c_T-1; t++){
			for(int i=0; i < c_dim; i++){
				REAL(tmp1)[0] = REAL(VECTOR_ELT(gamma,d))[i  + t*c_dim];
				REAL(tmp2)[0] = REAL(gamma_t)[i];
				REAL(gamma_t)[i] = REAL(_log_sum(tmp1,tmp2,base))[0];
				for(int j=0; j < c_dim; j++){
					REAL(tmp1)[0] = REAL(VECTOR_ELT(VECTOR_ELT(xi,d),t))[i + j*c_dim];
					REAL(tmp2)[0] = REAL(xi_t)[i + j*c_dim];
					REAL(xi_t)[i + j*c_dim] = REAL(_log_sum(tmp1,tmp2,base))[0];
				}
			}
		}
		SET_VECTOR_ELT(xi_sum_t,d,xi_t);
		SET_VECTOR_ELT(gamma_sum_t,d,gamma_t);
		UNPROTECT(2); /* xi_t, gamma_t */
	}
	
	/* combine results from different observation sequences */
	/* this uses equal weights for each sequence. There are other 
	 * (and possibly better) ways to do this */
	PROTECT(gamma_d = allocVector(REALSXP,c_dim));
	PROTECT(xi_d = allocMatrix(REALSXP,c_dim,c_dim));
	/* initialise */
	for(int i=0; i < c_dim; i++){
		REAL(gamma_d)[i] = REAL(VECTOR_ELT(gamma_sum_t,0))[i];
		REAL(pi)[i] = REAL(VECTOR_ELT(gamma,0))[i];
		for(int j=0; j < c_dim; j++){
			REAL(xi_d)[i + j*c_dim] = REAL(VECTOR_ELT(xi_sum_t,0))[i + j*c_dim];
		}	
	}
	/* sum over d */	
	for(int d=1; d < c_D; d++){
		for(int i=0; i < c_dim; i++){
			REAL(tmp1)[0] = REAL(VECTOR_ELT(gamma_sum_t,d))[i];
			REAL(tmp2)[0] = REAL(gamma_d)[i];
			REAL(gamma_d)[i] = REAL(_log_sum(tmp1,tmp2,base))[0];
			
			REAL(tmp1)[0] = REAL(VECTOR_ELT(gamma,d))[i];
			REAL(tmp2)[0] = REAL(pi)[i];
			REAL(pi)[i] = REAL(_log_sum(tmp1,tmp2,base))[0];
			
			for(int j=0; j < c_dim; j++){
				REAL(tmp1)[0] = REAL(VECTOR_ELT(xi_sum_t,d))[i + j*c_dim];
				REAL(tmp2)[0] = REAL(xi_d)[i + j*c_dim];
				REAL(xi_d)[i + j*c_dim] = REAL(_log_sum(tmp1,tmp2,base))[0];
			}	
		}
	}
	
	/* estimate new parameter values */
	PROTECT(trans = allocMatrix(REALSXP,c_dim,c_dim));
	for(int i=0; i < c_dim; i++){
		REAL(pi)[i] = (exp(REAL(pi)[i]) + REAL(init_prior)[i])/(c_D + c_sum);
			
		for(int j=0; j < c_dim; j++){
			REAL(trans)[i + j*c_dim] = (exp(REAL(xi_d)[i + j*c_dim]) + REAL(trans_prior)[i + j*c_dim]) / (exp(REAL(gamma_d)[i]) + REAL(prior_sum)[i]);
			//REAL(trans)[i + j*c_dim] = exp(REAL(trans)[i + j*c_dim]);	
		}
	}

	/* create list for return values */
	PROTECT(ret = allocVector(VECSXP,3));
	PROTECT(ret_names = allocVector(STRSXP,3));
	SET_STRING_ELT(ret_names,0,mkChar("pi"));
	SET_STRING_ELT(ret_names,1,mkChar("transition"));
	SET_STRING_ELT(ret_names,2,mkChar("gamma"));
	setAttrib(ret,R_NamesSymbol, ret_names);
	
	SET_VECTOR_ELT(ret,0,pi);
	SET_VECTOR_ELT(ret,1,trans);
	SET_VECTOR_ELT(ret,2,gamma);
	
	UNPROTECT(15); /* base, xi, log_trans, gamma, xi_sum_t, gamma_sum_t, tmp1, 
					  tmp2, gamma_d, xi_d, pi, trans, ret, ret_names, prior_sum 
					  (obs_prob) */

	return ret;
}

/* Viterbi algorithm */
SEXP _viterbi(SEXP hmm, SEXP obs){
	int c_states, c_obs_len = length(obs);
	SEXP dpm, trace, trans, trans_slot, obs_prob, init, init_slot, this_col, prob; 
	SEXP state_seq,ret, ret_names, hmm_trans;
	double *max_entry, c_init;

	/* get log of transition probabilities */
	PROTECT(trans_slot = allocVector(STRSXP,1));
	SET_STRING_ELT(trans_slot,0,mkChar("transition.matrix"));
	hmm_trans = R_do_slot(hmm,trans_slot);
	UNPROTECT(1);   /* trans_slot*/
	c_states = (int)sqrt(length(hmm_trans));
	PROTECT(trans = allocMatrix(REALSXP,c_states,c_states));
	for(int i=0; i < c_states * c_states; i++){
		REAL(trans)[i] = REAL(hmm_trans)[i];
		REAL(trans)[i] = log(REAL(trans)[i]);
	}

	
	/* calculate probability of most likely state sequence, given obs */
	PROTECT(dpm = allocMatrix(REALSXP,c_states, c_obs_len));
	PROTECT(trace = allocMatrix(INTSXP,c_states, c_obs_len-1));
	PROTECT(this_col = allocVector(REALSXP,c_states));
	
	/* initialisation */
	PROTECT(obs_prob = __get_obs_prob(hmm,REAL(obs)[0]));
	PROTECT(init_slot = allocVector(STRSXP,1));
	SET_STRING_ELT(init_slot,0,mkChar("init"));
	PROTECT(init = R_do_slot(hmm,init_slot));
	for(int i=0; i < c_states; i++){
		c_init = REAL(init)[i];
		REAL(dpm)[i] = log(c_init) + REAL(obs_prob)[i];
	}
	UNPROTECT(3);   /* init, obs_prob, init_slot */
	
	/* induction */
	for(int t=1; t < c_obs_len; t++){
		PROTECT(obs_prob = __get_obs_prob(hmm,REAL(obs)[t]));
		for(int i=0; i < c_states; i++){
			for(int j=0; j < c_states; j++){
				REAL(this_col)[j] = REAL(dpm)[j+(t-1)*c_states] + REAL(trans)[j+i*c_states];
			}
			max_entry = max(REAL(this_col),c_states);
			REAL(dpm)[i + t*c_states] = max_entry[0] + REAL(obs_prob)[i];
			INTEGER(trace)[i + (t-1)*c_states] = (int)max_entry[1];
			free(max_entry);
		}
		UNPROTECT(1);  /* obs_prob */
	}
	
	/* get the probability */
	PROTECT(prob = allocVector(REALSXP,1));
	max_entry = max(REAL(dpm)+c_states*(c_obs_len-1),c_states);
	REAL(prob)[0] = max_entry[0];
	
	/* calculate most likely state sequence given obs */
	PROTECT(state_seq = allocVector(INTSXP,c_obs_len));
	INTEGER(state_seq)[c_obs_len-1] = (int)max_entry[1]+1;
	free(max_entry);
	for(int t=c_obs_len-2; t >= 0; t--){
		INTEGER(state_seq)[t] = INTEGER(trace)[INTEGER(state_seq)[t+1]-1 + t*c_states]+1;	
	}
	
	/* create list for return values */
	PROTECT(ret = allocVector(VECSXP,3));
	PROTECT(ret_names = allocVector(STRSXP,3));
	SET_STRING_ELT(ret_names,0,mkChar("stateSeq"));
	SET_STRING_ELT(ret_names,1,mkChar("logProb"));
	SET_STRING_ELT(ret_names,2,mkChar("matrix"));
	setAttrib(ret,R_NamesSymbol, ret_names);
	
	SET_VECTOR_ELT(ret,0,state_seq);
	SET_VECTOR_ELT(ret,1,prob);
	SET_VECTOR_ELT(ret,2,dpm);	
	
	UNPROTECT(8);   /* trans, dpm, trace, this_col, prob, state_seq, ret, ret_names */
	
	return ret;
}

/* compute tau for t distribution EM */
SEXP __calcTau(SEXP gamma, SEXP obs, SEXP hmm){
	int c_states, c_obs_len = length(obs);
	SEXP tau, obs_prob, tau_sum, trans_slot, hmm_trans, base, tmp_obs_prob, tmp_tau_sum;
	
	/* get number of states */
	PROTECT(trans_slot = allocVector(STRSXP,1));
	SET_STRING_ELT(trans_slot,0,mkChar("transition.matrix"));
	hmm_trans = R_do_slot(hmm,trans_slot);
	UNPROTECT(1);   /* trans_slot*/
	c_states = (int)sqrt(length(hmm_trans));
	
	PROTECT(tau = allocMatrix(REALSXP,c_states,c_obs_len));
	PROTECT(tau_sum = allocVector(REALSXP, c_obs_len));
	PROTECT(base = allocVector(INTSXP,1));
	INTEGER(base)[0] = 0;
	PROTECT(tmp_obs_prob = allocVector(REALSXP,1));
	PROTECT(tmp_tau_sum = allocVector(REALSXP,1));
	
	/* calculate tau */
	for(int t=0; t < c_obs_len; t++){
		PROTECT(obs_prob = __get_obs_prob(hmm, REAL(obs)[t]));
		REAL(tau_sum)[t] = -INFINITY;
		/* get sum accross states for denominator */
		for(int i=0; i < c_states; i++){
			REAL(tau)[i + t*c_states] = REAL(obs_prob)[i] + REAL(gamma)[i + t*c_states];
			REAL(tmp_tau_sum)[0] = REAL(tau_sum)[t];
			REAL(tmp_obs_prob)[0] = REAL(tau)[i + t*c_states];
			REAL(tau_sum)[t] = REAL(_log_sum(tmp_tau_sum, tmp_obs_prob, base))[0];
		}
		/* combine everything to get tau */
		for(int i=0; i < c_states; i++){
			REAL(tau)[i + t*c_states] = REAL(tau)[i + t*c_states] - REAL(tau_sum)[t];	
		}
		UNPROTECT(1); /* obs_prob */
	}
	
	UNPROTECT(5); /* tau, tau_sum, base, tmp_obs_prob, tmp_tau_sum */
	return tau;		
}

/* Calculates probability that observation was produced by each state of the HMM.
 * Only t distributions are supported at the moment.*/
SEXP __get_obs_prob(SEXP hmm, double obs){
	double c_df,c_ncp,c_scale;
	int c_state_cnt;
	SEXP emission, prob, em_slot, comp_slot, em_class;
	
	PROTECT(em_slot = allocVector(STRSXP,1));
	SET_STRING_ELT(em_slot,0,mkChar("emission"));
	PROTECT(emission = R_do_slot(hmm,em_slot));
	c_state_cnt = length(emission);
	PROTECT(prob = allocVector(REALSXP,c_state_cnt));
	
	PROTECT(comp_slot = allocVector(STRSXP,1));
	SET_STRING_ELT(comp_slot,0,mkChar("components"));
	
	PROTECT(em_class = getAttrib(emission, R_ClassSymbol));
	
	for(int i=0;i<c_state_cnt;i++){
		/* Get mean and degrees of freedom for each component */
		c_df = REAL(R_do_slot(VECTOR_ELT(emission,i),comp_slot))[3];
		c_ncp = REAL(R_do_slot(VECTOR_ELT(emission,i),comp_slot))[1];
		c_scale = REAL(R_do_slot(VECTOR_ELT(emission,i),comp_slot))[2];
		REAL(prob)[i] = dt((obs - c_ncp)/sqrt(c_scale),c_df,1) - 0.5*log(c_scale);	
	}
	
	UNPROTECT(5);
	return prob;
}
SEXP _get_obs_prob(SEXP hmm, SEXP obs){
	return __get_obs_prob(hmm,REAL(obs)[0]);
}
