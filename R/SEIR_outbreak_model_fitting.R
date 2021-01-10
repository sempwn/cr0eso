
#' list of parameterisations for priors
#' @export
prior_list <- list(gamma_mean = 0.125,gamma_sd = 0.0125,
                   sigma_mean = 0.2,sigma_sd = 0.025,
                   S0_mean = 0.9, S0_sd = 0.01,
                   r0_mean=3.0,r0_sd=1.0,
                   zeta_mean=0.1,zeta_sd=0.1,
                   tau_mean = 0.1)


#' SEIR model fitting method to a multiple outbreak model
#' @description
#' Create an instance of the hierarchical SEIR stan model incorporating
#' various data elements and sample model.
#' @param tmax Total number of time-points in observation
#' @param n_outbreaks Total number of outbreaks
#' @param outbreak_cases Numnber of daily reported cases by outbreak
#' @param outbreak_sizes  The total size of each facility (initial number of suscepitble and exposed)
#' @param intervention_switch Describes whether interventions occur in data (default TRUE)
#' @param multilevel_intervention Describes whether intervention occurs
#' @param prior_list List of priors. See `prior_list`
#' @param iter number of iterations of MCMC
#'
#' @examples
#' tmax <- 5
#' pop_size <- 100
#' dim(pop_size) <- c(1)
#' example_incidence <- matrix(c(1,1,2,3,2),ncol=1)
#' fit <- seir_model_fit(tmax,1,example_incidence,pop_size)
#' @author Mike Irvine
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' @export
seir_model_fit <- function(tmax,
                           n_outbreaks,
                           outbreak_cases,
                           outbreak_sizes,
                           intervention_switch = TRUE,
                           multilevel_intervention = FALSE,
                           priors = prior_list,
                           iter=600
                           ){

  # define STAN data
  stan_data = list(n_obs = tmax,
                   n_outbreaks = n_outbreaks,
                   n_difeq = 5,
                   n_fake = tmax,
                   y = outbreak_cases,
                   intervention_switch = TRUE,
                   multilevel_intervention = FALSE,
                   independent_r0 = FALSE,
                   independent_zeta = FALSE,
                   n_prior_mean = outbreak_sizes,
                   tau_prior_mean = priors$tau_mean,
                   t0 = 0,
                   tn = tmax,
                   ts = c(1:tmax),
                   fake_ts = c(1:tmax),
                   # define priors
                   gamma_mean = priors$gamma_mean,
                   gamma_sd = priors$gamma_sd,
                   sigma_mean = priors$sigma_mean,
                   sigma_sd = priors$sigma_sd,
                   S0_mean = priors$S0_mean,
                   S0_sd = priors$S0_sd,
                   r0_mean = priors$r0_mean,
                   r0_sd = priors$r0_sd,
                   zeta_mean = priors$zeta_mean,
                   zeta_sd = priors$zeta_sd
  )

  # Which parameters to monitor in the model:
  params_monitor = c("y0", "params","r0k" , "r0" , "r0_sigma", "zetak",
                     "incidence" , "hyper_priors","fake_incidence",
                     "predictive_r0")


  # Fit and sample from the posterior
  mod = rstan::sampling(stanmodels$hierarchical_SEIR_incidence_model,
             data = stan_data,
             pars = params_monitor,
             seed = 42,
             chains = 4,
             warmup = floor(iter/2),
             iter = iter,
             control = list(adapt_delta = 0.95))

  list(model = mod)
}
