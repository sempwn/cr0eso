
#' list of parameterisations for priors
#' @export
prior_list <- list(gamma_mean = 0.125,gamma_sd = 0.0125,
                   sigma_mean = 0.2,sigma_sd = 0.025,
                   S0_mean = 0.9, S0_sd = 0.01,
                   r0_mean=3.0,r0_sd=1.0,
                   zeta_mean=0.1,zeta_sd=0.1,
                   tau_mean = 0.1,
                   phi = c(1, 1))


#' SEIR model fitting method to a multiple outbreak model
#' @description
#' Create an instance of the hierarchical SEIR Stan model incorporating
#' various data elements and sample model.
#' @param stan_model [rstan] model object
#' @param tmax Total number of time-points in observation
#' @param n_outbreaks Total number of outbreaks
#' @param outbreak_cases Number of daily reported cases by outbreak
#' @param outbreak_sizes  The total size of each facility (initial number of susceptible and exposed)
#' @param intervention_switch Describes whether interventions occur in data (default TRUE)
#' @param multilevel_intervention Describes whether intervention occurs
#' @param prior_list List of priors. See [prior_list]
#' @param chains Number of chains to sample
#' @param iter number of iterations of MCMC
#' @param seed The seed for random number generation. Set to replicate results.
#' @param fit_type string "NUTS" or "VB" (VB quicker but less accurate).
#' @param data_model string "poisson" or "negative_binomial". If negative_binomial
#'   selected then uses phi prior to control for overdispersion
#'
#' @examples
#' stan_mod <- rstan::stan_model(system.file("stan",
#'   "hierarchical_SEIR_incidence_model.stan", package = "cr0eso"))
#' tmax <- 5
#' pop_size <- 100
#' dim(pop_size) <- c(1)
#' example_incidence <- matrix(c(1,1,2,3,2),ncol=1)
#' fit <- seir_model_fit(stan_model = stan_mod, tmax,1,example_incidence,pop_size)
#' @author Mike Irvine
#' @return An object of class `stanfit` returned by \link[rstan]{sampling}
#' @export
seir_model_fit <- function(
                           stan_model = NULL,
                           tmax,
                           n_outbreaks,
                           outbreak_cases,
                           outbreak_sizes,
                           intervention_switch = TRUE,
                           multilevel_intervention = FALSE,
                           priors = prior_list,
                           chains = 4,
                           iter=600,
                           seed = 42,
                           fit_type = "NUTS",
                           data_model="poisson"
                           ){

  if(is.null(stan_model)){
    stop("Need to include a stan model. Check examples")
  }

  # parameter checks
  if(!fit_type %in% c("NUTS","VB")){
    stop(glue::glue("{fit_type} should be NUTS or VB."))
  }

  if(!data_model %in% c("poisson","negative_binomial")){
    stop(glue::glue("{data_model} should be poisson or negative_binomial."))
  }

  data_model_code <- dplyr::case_when(data_model == "poisson" ~ 0,
                                      data_model == "negative_binomial" ~ 1)
  # define STAN data
  stan_data = list(n_obs = tmax,
                   n_outbreaks = n_outbreaks,
                   n_difeq = 5,
                   n_fake = tmax,
                   y = outbreak_cases,
                   intervention_switch = intervention_switch,
                   multilevel_intervention = multilevel_intervention,
                   independent_r0 = FALSE,
                   independent_zeta = FALSE,
                   data_model = data_model_code,
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
                   zeta_sd = priors$zeta_sd,
                   phi_prior = priors$phi
  )

  # Which parameters to monitor in the model:
  params_monitor = c("y0", "params","r0k" , "r0" , "r0_sigma", "zetak",
                     "incidence" , "hyper_priors","counterfactual_cases",
                     "pp_cases","predictive_r0","phi")


  # Fit and sample from the posterior
  if(fit_type=="NUTS"){
    mod <- rstan::sampling(stan_model,
               data = stan_data,
               pars = params_monitor,
               seed = seed,
               chains = chains,
               warmup = floor(iter/2),
               iter = iter,
               control = list(adapt_delta = 0.95))
  }else{

      cat("Sampling with the VB algorithm.\n")
      mod <- rstan::vb(
        stanmodels$hierarchical_SEIR_incidence_model,
        data = stan_data,
        iter = iter,
        seed = seed,
        pars = params_monitor,
        algorithm = "fullrank",
        tol_rel_obj = 0.0001,
        importance_resampling = TRUE
      )

  }

  list(model = mod)
}
