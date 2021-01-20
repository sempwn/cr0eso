// STAN hierarchical model of SEIR outbreaks
//


functions {

  real[] SEI(real t, // time
            real[] y, // state
            real[] params, // parameters (involved in fitting)
            real[] x_r, // data (real)
            int[] x_i) { // data (integer)
      //define states
      real S = y[1];
      real E = y[2];
      real I = y[3];
      real R = y[4];

      // define params
      real r0 = params[1];
      real sigma = params[2];
      real gamma = params[3];
      real n = params[4];
      real tau = params[5];
      real zeta = params[6];

      int intervention_switch = x_i[1];

      //define placeholder values
      real r;

      real dydt[5];

      // modified r0 dependent on intervention
      if(intervention_switch != 0){
        r = r0*(t<=tau) +  (t>tau)*r0*exp(zeta*(tau-t));
      }else{
        r = r0;
      }

      dydt[1] = - r * gamma * S * I; //dS
      dydt[2] = r * gamma * S * I - sigma * E; //dE
      dydt[3] = sigma * E - gamma * I; //dI
      dydt[4] = gamma * I; //dR
      dydt[5] = n * sigma * E; //d incidence

      return dydt;
    }


    // simulate random numbers from a truncated normal distribution
    // for drawing from priors
    real normal_lb_rng(real mu, real sigma, real lb) {
        real p = normal_cdf(lb, mu, sigma);  // cdf for bounds
        real u = uniform_rng(p, 1);
        return (sigma * inv_Phi(u)) + mu;  // inverse cdf for value
    }

}

data {
  int<lower = 1> n_obs; // Number of days sampled
  int<lower = 1> n_outbreaks; // Number of outbreaks observed
  int<lower = 1> n_difeq; // Number of differential equations in the system
  int<lower = 1> n_fake; // This is to generate "predicted"/"unsampled" data
  int<lower=0, upper=1> intervention_switch; // logical: include effects of intervention
  int<lower=0, upper=1> multilevel_intervention; //logical: if included intervention strength zeta is modelled as a
                                                 //         hierarchical parameter where each location has an individual zeta
  int<lower=0, upper=1> independent_r0; //logical: if true then each R0 modelled independently
                                        // useful for creating shrinkage plots

  int<lower=0, upper=1> independent_zeta; //logical: if true then each R0 modelled independently
                                        // useful for creating shrinkage plots

  int y[n_obs,n_outbreaks]; // The binomially distributed data

  // priors
  real<lower=0> n_prior_mean[n_outbreaks]; // mean of population size
  real<lower=0> tau_prior_mean; // mean of tau (start of intervention)
  real<lower=0> gamma_mean;
  real<lower=0> gamma_sd;
  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;
  real<lower=0> S0_mean;
  real<lower=0> S0_sd;
  real<lower=0> r0_mean;
  real<lower=0> r0_sd;
  real<lower=0> zeta_mean;
  real<lower=0> zeta_sd;


  //time-points
  real t0; // Initial time point (zero)
  real tn; // final time-point
  real ts[n_obs]; // Time points that were sampled

  real fake_ts[n_fake]; // Time points for "predicted"/"unsampled" data
}

transformed data {
  real x_r[0];
  int x_i[1];

  x_i[1] = intervention_switch;
}

parameters {


  // priors (will be updated)
  real<lower=0> r0;
  real<lower=0> r0_sigma;

  // multilevel r0 (non-centred parameterization)
  vector<lower=0>[n_outbreaks] r0k_raw;

  real std_gamma;
  real std_sigma;
  real<lower=0> tau;
  real<lower=0> zeta;

  // multilevel intervention parameters
  real<lower=0> zeta_sigma;
  vector<lower=0>[n_outbreaks] zetak_raw;

  real<lower = 0> n[n_outbreaks]; //effective population size
  real std_S0; // Initial fraction of hosts susceptible

}

transformed parameters{
  real y_hat[n_obs, n_difeq]; // Output from the ODE solver
  real incidence[n_obs,n_outbreaks];
  real y0[n_difeq]; // Initial conditions for both S and I
  vector<lower=0>[n_outbreaks] r0k;
  vector<lower=0>[n_outbreaks] zetak;
  real params[6]; // Model parameters
  real<lower=0> gamma; //rate (days^-1) version of inv_gamma (recovery period in days)
  real<lower=0> sigma; //rate version of inv_sigma (incubation period in days)

  real<lower = 0, upper = 1> S0; // Initial fraction of hosts susceptible

  real eps = 1e-4;

  // se rates from their periods
  gamma = gamma_mean + gamma_sd*std_gamma;
  sigma = sigma_mean + sigma_sd*std_sigma;

  // S0 from standardized S0
  S0 = S0_mean + S0_sd*std_S0;

  // define noncentred r0k
  if(independent_r0==1){
    r0k = r0_mean + r0_sd * r0k_raw;
  }else{
    r0k = r0 + r0_sigma * r0k_raw;
  }

  //define noncentred zetak
  if(independent_zeta == 1){
    zetak = zeta_mean + zeta_sd * zetak_raw;
  } else{
    zetak = zeta + zeta_sigma * zetak_raw;
  }


  for(k in 1:n_outbreaks){
    // create posterior vector
    // multilevel parameter
    params[1] = r0k[k];
    params[2] = sigma;
    params[3] = gamma;
    params[4] = n[k];

    // if including intervention has int params
    // else hard set to zero so does not effect R0
    if(intervention_switch != 0){
      params[5] = tau;
      // if no multilevel intervention draw from general distribution
      // otherwise draw from individual zetak distribution
      if(multilevel_intervention == 0){ // SA: may create sampling problems since the other still exists
        params[6] = zeta;
      }else{
        params[6] = zetak[k];
      }

    }else{ // SA: This may create sampling problems since zeta and zetak still exist:
      params[5] = 0;
      params[6] = 0;
    }



    y0[1] = S0; //S
    y0[2] = 1 - S0; //E
    y0[3] = 0; //I
    y0[4] = 0; //R
    y0[5] = 0; //cumulative incidence

    y_hat = integrate_ode_rk45(SEI, y0, t0, ts, params, x_r, x_i);
    incidence[,k] = y_hat[,5];
  }

}

model {

  // list priors that will be updated
  r0 ~ normal(r0_mean, r0_sd); // multilevel prior previous var 1.0
  r0_sigma ~ normal(0,1);
  r0k_raw ~ std_normal();

  std_gamma ~ std_normal();
  std_sigma ~ std_normal();
  n ~ normal(n_prior_mean,10);

  std_S0 ~ std_normal();

  // if including intervention has int params
  // else hard set to zero so does not effect R0
  tau ~ normal(tau_prior_mean,1.0);
  zeta ~ normal(zeta_mean,zeta_sd); // previous var 0.1
  // multilevel priors for zeta. Only used if multilevel priors are switched on
  zeta_sigma ~ normal(0,1);
  zetak_raw ~ std_normal();
  //zetak ~ normal(zeta, zeta_sigma);
  for(k in 1:n_outbreaks){
    for(i in 2:n_obs){
      y[i,k] ~ poisson(max([incidence[i,k] - incidence[i-1,k], 1e-10 ]));
    }
  }


}

generated quantities {

  // Generate predicted data over the whole time series:
  real fake_I[n_fake, n_difeq];

  // generate posterior predictive check for cases
  real pp_cases[n_obs,n_outbreaks];

  // Generate posterior predictive check for counterfactual
  real cparams[6]; // Model parameters
  real counterfactual_cases[n_obs,n_outbreaks];
  real diff_I;

  // Model prior parameters
  real hyper_priors[5];
  real<lower=0> p_r0;
  real<lower=0> p_gamma;
  real<lower=0> p_sigma;
  real<lower=0> p_tau;
  real<lower=0> p_zeta;
  //R0 predictive distribution
  real<lower=0> predictive_r0;


  // posterior predictive cases counterfactual scenario
  for(k in 1:n_outbreaks){
    // create posterior vector
    // multilevel parameter
    cparams[1] = r0k[k];
    cparams[2] = sigma;
    cparams[3] = gamma;
    cparams[4] = n[k];

    // set intervention effect to 0
    cparams[5] = 0;
    cparams[6] = 0;

    fake_I = integrate_ode_rk45(SEI, y0, t0, ts, cparams, x_r, x_i);

    counterfactual_cases[1,k] = 0;
    for(i in 2:n_obs){
      diff_I = fake_I[i,5] - fake_I[i-1,5]; //y_hat[,5] cumulative incidence
      if(diff_I < 0){
        diff_I = 1e-3; // small value if incidence is negative.
      }
      counterfactual_cases[i,k] = poisson_rng(diff_I);
    }
  }



  // posterior predictive cases
  for(k in 1:n_outbreaks){
    // set day 0 cases to 0
    pp_cases[1,k] = 0;
    for(i in 2:n_obs){
      pp_cases[i,k] = poisson_rng(max([incidence[i,k] - incidence[i-1,k], 1e-10 ]));
    }
  }

  // generate predictive distribution for R0
  predictive_r0 = normal_lb_rng(r0,r0_sigma,0);

  // generate priors for testing purposes

  //priors
  p_r0 = normal_lb_rng(r0_mean, r0_sd, 0);
  p_gamma = normal_lb_rng(gamma_mean,gamma_sd, 0);
  p_sigma = normal_lb_rng(sigma_mean,sigma_sd, 0);
  p_tau = normal_lb_rng(tau_prior_mean,1.0, 0);
  p_zeta = normal_lb_rng(zeta_mean,zeta_sd, 0);


  // create vector of priors
  hyper_priors[1] = p_r0;
  hyper_priors[2] = p_sigma;
  hyper_priors[3] = p_gamma;
  hyper_priors[4] = p_tau;
  hyper_priors[5] = p_zeta;


}





