// Model asocial reinforcement learning 
// (Rescorla-Wagner updating + Softmax choice)
// Assumtion: alpha and beta may be correlated

data {
  int<lower = 0> All;  // number of observations
  int<lower = 0> Nsub;  // number of subjects
  int<lower = 0> Ncue;  // number of options
  int<lower = 0> Ntrial;  // number of trials per subject
  int<lower = 0> sub[All];  // subject index
  int<lower = 0> Y[All];  // index of chosen option: 0 => missing ()
  int<lower = 0> trial[All];  // trial number
  real payoff[All];  // payoff
}

parameters {
  //vector[Nsub] alpha_raw;  // learning rate -- raw valuable
  //vector[Nsub] beta_raw;  // softmax parameter
  real mu_alpha;
  real mu_beta;
  real<lower=0> s_alpha;
  real<lower=0> s_beta;
  
  // Varying effects clustered on individual, defined as deviations from grand mean
  // We do the non-centered version with the Cholesky factorization

  matrix[2,Nsub] z_ID;               //Matrix of uncorrelated z - values
  vector<lower=0>[2] sigma_ID;       //SD of parameters among individuals
  cholesky_factor_corr[2] L_Rho_ID;    // This is the Cholesky factor: if you multiply this matrix and it's transpose you get correlation matrix

}

transformed parameters {
  matrix[Nsub,2] v_ID; // Matrix of varying effects for each individual for alpha, beta
  
  matrix[Ncue, Ntrial] Q[Nsub];  // Q-values for each target
  matrix[Ncue, Ntrial] q[Nsub]; // softmax choice (log-)probability transformed from Q values
  vector[Nsub] counter;
  vector[Nsub] logit_alpha;  // learning rate -- logit scale
  vector[Nsub] log_beta;  // softmax parameter -- log scale
  vector<lower = 0, upper = 1>[Nsub] alpha;  // learning rate
  vector<lower = 0>[Nsub] beta;  // inverse temperature (softmax)
  
  v_ID = ( diag_pre_multiply( sigma_ID , L_Rho_ID ) * z_ID )'; // the ' on the end of this line is a transposing function

  counter = rep_vector(0, Nsub); // tracking each subject's first experience timing

  for(i in 1:Nsub) {
    logit_alpha[i] = mu_alpha + s_alpha * v_ID[i, 1];
    log_beta[i] = mu_beta + s_beta * v_ID[i, 2];
    //logit_alpha[i] = mu_alpha + s_alpha * alpha_raw[i];
    //log_beta[i] = mu_beta + s_beta * beta_raw[i];
    beta[i] = exp(log_beta[i]);
    alpha[i] = 1/(1+exp(-logit_alpha[i]));
  }

  for(idx in 1:All) {
    // SETTING INITIAL Q-VALUE
    if(trial[idx] == 1)
      Q[sub[idx]][1:Ncue, trial[idx]] = rep_vector(0, Ncue);
    // SETTING INITIAL Q-VALUE -- END

    // CHOICE PROBABILITY
    q[sub[idx]][1:Ncue, trial[idx]] = 
      log_softmax( 
        Q[sub[idx]][, trial[idx]] * beta[sub[idx]] 
        );
    
    //for(c in 1:Ncue) {
        // Choice probability by Asocial softmax rule
      //q[sub[idx]][c, trial[idx]] = log_softmax( Q[sub[idx]][, trial[idx]] * beta[sub[idx]] )[c];
    //}
    
    // CHOICE PROBABILITY -- END

    // Q-VALUE UPDATE
    if(trial[idx] < Ntrial) {
      
      Q[sub[idx]][, (trial[idx]+1)] = Q[sub[idx]][, trial[idx]];
      //for(c in 1:Ncue)
        //Q[sub[idx]][c, (trial[idx]+1)] = Q[sub[idx]][c, trial[idx]];
        
      // update chosen option
      if(Y[idx]>0) {
        if(counter[sub[idx]]==0)
          { // Updating all Q-values at the 1st experience
            Q[sub[idx]][, (trial[idx]+1)] = (1-alpha[sub[idx]])*Q[sub[idx]][, trial[idx]] + alpha[sub[idx]]*payoff[idx];
            //for(c in 1:Ncue)
              //Q[sub[idx]][c, (trial[idx]+1)] = (1-alpha[sub[idx]])*Q[sub[idx]][c, trial[idx]] + alpha[sub[idx]]*payoff[idx];
            counter[sub[idx]] = 1;
          }
        else
          { // Updating chosen option's Q-value
            Q[sub[idx]][Y[idx], (trial[idx]+1)] =
              (1-alpha[sub[idx]])*Q[sub[idx]][Y[idx], trial[idx]] + alpha[sub[idx]]*payoff[idx];
          }
      }
    }
    // Q-VALUE UPDATE -- END
  }
}

model {
  // position and scales of the learing parameters
  mu_alpha ~ std_normal(); // prior
  mu_beta ~ normal(-1, 1); // prior
  s_alpha ~ exponential(1); //std_normal(); // prior
  s_beta ~ exponential(1); //std_normal(); // prior
  //alpha_raw ~ student_t(4, 0, 1); // prior
  //beta_raw ~ student_t(4, 0, 1); // prior
  
  // Varying effects
  to_vector(z_ID) ~ std_normal();
  sigma_ID ~ exponential(1);
  L_Rho_ID ~ lkj_corr_cholesky(2); // A weakly informative prior for the correlation matrix

  for(idx in 1:All) {
    if(Y[idx] > 0) {
      target += categorical_lpmf( Y[idx] | exp(q[sub[idx]][,trial[idx]]) );
    }
  }
}

generated quantities {
  vector[Nsub] log_lik;
  matrix[2,2] Rho_ID; // the correlation matrix
  
  Rho_ID = multiply_lower_tri_self_transpose(L_Rho_ID);
  
  log_lik = rep_vector(0, Nsub); // initial values for log_lik
  for(idx in 1:All) {
    if( Y[idx] > 0 ) {
      log_lik[sub[idx]] = log_lik[sub[idx]] + q[sub[idx]][Y[idx],trial[idx]];
    }
  }
}

