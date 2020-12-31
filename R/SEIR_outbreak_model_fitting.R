

#' SEIR model fitting method to a multiple outbreak model
#' @description
#' Create an instance of the hierarchical SEIR stan model incorporating
#' various data elements and sample model.
#' @param tmax - total number of time-points in observation
#' @param n_outbreaks - total number of outbreaks
#' @param outbreak_cases - Numnber of daily reported cases by outbreak
#' @param outbreak_sizes - The total size of each facility (initial number of suscepitble and exposed)
#' @param intervention_switch - Describes whether interventions occur in data (default TRUE)
#' @param multilevel_intervention - describes whether intervention occurs
#' @param iter number of iterations of MCMC
#'
#' @examples
#' tmax <- 30
#' example_incidence <- c(1,1,2,3,2)
#' seir_model_fit(tmax,c(1),example_incidence,c(100))
#' @author Mike Irvine
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' @export
seir_model_fit <- function(tmax,
                           n_outbreaks,
                           outbreak_cases,
                           outbreak_sizes,
                           intervention_switch = TRUE,
                           multilevel_intervention = FALSE,
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
                   tau_prior_mean = 0.1,
                   t0 = 0,
                   tn = tmax,
                   ts = c(1:tmax),
                   fake_ts = c(1:tmax)
  )

  # Which parameters to monitor in the model:
  params_monitor = c("y0", "params","r0k" , "r0" , "r0_sigma",
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
