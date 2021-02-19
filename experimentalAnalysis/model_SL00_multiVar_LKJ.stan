// SL00: without temporal dynamics in soc_weight and temperature
// Baseline AL: RW rule + Softmax choice
// 14/12/2020
// Wataru Toyokawa (wataru.toyokawa@uni-konstanz.de)
// Assumption: all learning parameters could be correlated each other.
// The non-centered variance is obtained from

functions {
  real net_choice_prob_lpmf(int Y, real soc, int trial, vector q, vector q_f, vector f) {
    if(trial == 1 || sum(f)==3e-10) {
      // At the first trial or when there is no one else, rely only on asocial learning
      return categorical_lpmf(Y | exp(q));
    }
    else {
      return log_sum_exp( log( 1-soc ) + categorical_lpmf(Y | exp(q)),
                          log( soc ) + categorical_lpmf(Y | exp(q_f))
                        );
    }
  }
}

data {
    int<lower = 0> All;  // number of observations
    int<lower = 0> Nsub;  // number of subjects
    int<lower = 0> Ncue;  // number of cues
    int<lower = 0> Ntrial;  // number of trials per subject
    int<lower = 0> sub[All];  // subject index
    int<lower = 0> Y[All];  // index of chosen option: 0 => missing
    int<lower = 0> trial[All];  // trial number
    real payoff[All];  // payoff
    real F[Nsub,Ncue,Ntrial];
}

transformed data {
  matrix[Ncue, Ntrial] f[Nsub]; // frequency information (transfered)
  vector<lower = 0>[All] trial_real;
  for(idx in 1:All) {
    trial_real[idx] = trial[idx];
    for(c in 1:Ncue) {
      f[sub[idx]][c,trial[idx]] = F[sub[idx]][c,trial[idx]] + 1e-01;
    }
  }
}

parameters {
    real mu_alpha;
    real mu_beta;
    real mu_soc0;
    real mu_theta;
    // real mu_soc_slope;
    // real mu_epsilon;

    real<lower=0> s_alpha;
    real<lower=0> s_beta;
    real<lower=0> s_theta;
    real<lower=0> s_soc0;
    // real<lower=0> s_soc_slope;
    // real<lower=0> s_epsilon;

    // Varying effects clustered on individual, defined as deviations from grand mean
    // We do the non-centered version with the Cholesky factorization

    matrix[4,Nsub] z_ID;               //Matrix of uncorrelated z - values
    vector<lower=0>[4] sigma_ID;       //SD of parameters among individuals
    cholesky_factor_corr[4] L_Rho_ID;    // This is the Cholesky factor:
    //if you multiply this matrix and it's transpose you get correlation matrix
}

transformed parameters {
  matrix[Nsub,4] v_ID; // Matrix of varying effects for each individual for alpha, beta

  matrix[Ncue, Ntrial] Q[Nsub];  // value function for each target
  matrix[Ncue, Ntrial] q[Nsub]; // softmax choice (log) probability
  matrix[Ncue, Ntrial] q_f[Nsub]; // frequency dependent copying (log) probability
  vector[Ntrial] Fsum[Nsub];
  simplex[Ncue] q_f_simplex[Nsub];

  vector[Nsub] logit_alpha;  // learning rate -- raw valuable
  vector[Nsub] log_beta;  // softmax parameter
  // vector[Nsub] epsilon;  // softmax parameter
  vector[Nsub] logit_soc0; // copying rate -- raw
  // vector[Nsub] logit_soc_slope;  // softmax parameter
  vector[Nsub] theta;  // conformity exponent

  vector<lower = 0, upper = 1>[Nsub] alpha;  // learning rate
  vector<lower = 0, upper = 1>[Nsub] soc;  // learning rate
  vector<lower = 0>[Nsub] beta;  // inverse temperature (softmax)
  vector[Nsub] counter;

  v_ID = ( diag_pre_multiply( sigma_ID , L_Rho_ID ) * z_ID )';
  // the ' on the end of this line is a transposing function

  counter = rep_vector(0, Nsub); // tracking each subject's first experience timing

  for(i in 1:Nsub) {
    logit_alpha[i] = mu_alpha + s_alpha * v_ID[i, 1];
    log_beta[i] = mu_beta + s_beta * v_ID[i, 2];
    // epsilon[i] = mu_epsilon + s_epsilon * v_ID[i, 3];
    logit_soc0[i] = mu_soc0 + s_soc0 * v_ID[i, 3]; // logit_soc0 is translated into soc in the final function
    // logit_soc_slope[i] = mu_soc_slope + s_soc_slope * v_ID[i, 5];
    theta[i] = mu_theta + s_theta * v_ID[i, 4];

    beta[i] = exp(log_beta[i]);
    alpha[i] = 1/(1+exp(-logit_alpha[i]));
    soc[i] = 1/(1+exp(-logit_soc0[i]));
  }

  for(idx in 1:All) {
    // === SETTING INITIAL Q-VALUE
    if(trial[idx] == 1)
      Q[sub[idx]][1:Ncue, trial[idx]] = rep_vector(0, Ncue);
    // === SETTING INITIAL Q-VALUE -- END

    // === DENOMINATOR OF SOCIAL FREQUENCY INFORMATION
    Fsum[sub[idx]][trial[idx]] = 0;
    for(c in 1:Ncue) {
      // f1^theta + f2^theta
      Fsum[sub[idx]][trial[idx]] = Fsum[sub[idx]][trial[idx]] + f[sub[idx]][c,trial[idx]]^theta[sub[idx]];
    }
    // === DENOMINATOR OF SOCIAL FREQUENCY INFORMATION -- END

    // === ASOCIAL AND SOCIAL CHOICE PROBABILITY
    q[sub[idx]][1:Ncue, trial[idx]] =
      log_softmax(
        Q[sub[idx]][, trial[idx]] * beta[sub[idx]]
        );
    q_f_simplex[sub[idx]] = exp(theta[sub[idx]]*log(f[sub[idx]][,trial[idx]])-log(Fsum[sub[idx]][trial[idx]]));
    q_f[sub[idx]][1:Ncue, trial[idx]] = log(q_f_simplex[sub[idx]]); // log transform
    //for(c in 1:Ncue) {
      // Choice probability by Asocial sofimax rule
      //q[sub[idx]][c, trial[idx]] =
        //log_softmax(
          //Q[sub[idx]][, trial[idx]] * exp( log_beta[sub[idx]] + epsilon[sub[idx]]*(trial_real[idx]/70) )
        //)[c];
      // social frequency influence
      //q_f_simplex[sub[idx]][c] = exp(theta[sub[idx]]*log(f[sub[idx]][c,trial[idx]])-log(Fsum[sub[idx]][trial[idx]]));
      //q_f[sub[idx]][c, trial[idx]] = log(q_f_simplex[sub[idx]][c]); // log transform
    //}
    // === ASOCIAL AND SOCIAL CHOICE PROBABILITY -- END

    // === Q-VALUE UPDATE
    if(trial[idx] < Ntrial) {
      Q[sub[idx]][, (trial[idx]+1)] = Q[sub[idx]][, trial[idx]];
      // update chosen option
      if(Y[idx]>0) {
        if( counter[sub[idx]]==0 )
          { // Updating all Q-values at the 1st experience
            Q[sub[idx]][, (trial[idx]+1)] =
              (1-alpha[sub[idx]])*Q[sub[idx]][, trial[idx]] + alpha[sub[idx]]*payoff[idx];
            counter[sub[idx]] = 1;
          }
        else
          { // Updating chosen option's Q-value
            Q[sub[idx]][Y[idx], (trial[idx]+1)] =
              Q[sub[idx]][Y[idx], trial[idx]] + alpha[sub[idx]] * (payoff[idx] - Q[sub[idx]][Y[idx], trial[idx]]);
              //(1-alpha[sub[idx]])*Q[sub[idx]][Y[idx], trial[idx]] + alpha[sub[idx]]*payoff[idx];
          }
      }
    }
    // === Q-VALUE UPDATE -- END
  }
}

model {
  mu_alpha ~ std_normal(); // for [0,1] parameter
  mu_beta ~ normal(-1, 1); // prior 1 ~ 5
  mu_soc0 ~ std_normal(); // prior -3 ~ 3
  mu_theta ~ std_normal(); // prior -3 ~ 7
  // mu_soc_slope ~ std_normal(); //normal(-2, 3); // prior -4 ~ 2
  // mu_epsilon ~ std_normal(); //normal(3, 2);

  s_alpha ~ exponential(1); // prior
  s_beta ~ exponential(1); // prior
  s_soc0 ~ exponential(1); // prior
  s_theta ~ exponential(1); // prior
  // s_soc_slope ~ exponential(1); // prior
  // s_epsilon ~ exponential(1); // prior

  /// Varying effects
  to_vector(z_ID) ~ std_normal();
  sigma_ID ~ exponential(1);
  L_Rho_ID ~ lkj_corr_cholesky(4); // A weakly informative prior for the correlation matrix 0~1: almost flat, >2: weak correlation

  for(idx in 1:All) {
    if(Y[idx] > 0) {
      target +=
        net_choice_prob_lpmf(Y[idx] | soc[sub[idx]], trial[idx], q[sub[idx]][,trial[idx]], q_f[sub[idx]][,trial[idx]], f[sub[idx]][,trial[idx]] );
    }
  }
}

generated quantities {
  vector[Nsub] log_lik;
  // vector<lower = 0, upper = 1>[Ntrial] soc[Nsub];
  // vector[Ntrial] netBeta[Nsub];
  matrix[4,4] Rho_ID; // the correlation matrix

  Rho_ID = multiply_lower_tri_self_transpose(L_Rho_ID);

  log_lik = rep_vector(0, Nsub); // initial values for log_lik
  // for (i in 1:Nsub) {
  //   for (t in 1:Ntrial) {
  //     // soc[i][t] = 1/(1+exp(-(logit_soc0[i]+(logit_soc_slope[i]*t/70))));
  //     // netBeta[i][t] = exp( log_beta[i] + epsilon[i]*t/70 );
  //   }
  // }
  for(idx in 1:All) {
    if(Y[idx] > 0) {
      if(trial[idx] == 1 || sum(f[sub[idx]][,trial[idx]])==3e-10)
        log_lik[sub[idx]] = log_lik[sub[idx]] + q[sub[idx]][Y[idx],trial[idx]];
      else
        log_lik[sub[idx]] = log_lik[sub[idx]] +
          log( ( 1-soc[sub[idx]] )*exp(q[sub[idx]][Y[idx],trial[idx]])+
                soc[sub[idx]]*exp(q_f[sub[idx]][Y[idx],trial[idx]])
          );
    }
  }
}

