% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEIR_outbreak_model_fitting.R
\name{seir_model_fit}
\alias{seir_model_fit}
\title{SEIR model fitting method to a multiple outbreak model}
\usage{
seir_model_fit(
  stan_model = NULL,
  tmax,
  n_outbreaks,
  outbreak_cases,
  outbreak_sizes,
  intervention_switch = TRUE,
  multilevel_intervention = FALSE,
  priors = prior_list,
  chains = 4,
  iter = 600,
  seed = 42,
  fit_type = "NUTS",
  data_model = "poisson"
)
}
\arguments{
\item{stan_model}{[rstan] model object}

\item{tmax}{Total number of time-points in observation}

\item{n_outbreaks}{Total number of outbreaks}

\item{outbreak_cases}{Number of daily reported cases by outbreak}

\item{outbreak_sizes}{The total size of each facility (initial number of susceptible and exposed)}

\item{intervention_switch}{Describes whether interventions occur in data (default TRUE)}

\item{multilevel_intervention}{Describes whether intervention occurs}

\item{chains}{Number of chains to sample}

\item{iter}{number of iterations of MCMC}

\item{seed}{The seed for random number generation. Set to replicate results.}

\item{fit_type}{string "NUTS" or "VB" (VB quicker but less accurate).}

\item{data_model}{string "poisson" or "negative_binomial". If negative_binomial
selected then uses phi prior to control for overdispersion}

\item{prior_list}{List of priors. See [prior_list]}
}
\value{
An object of class `stanfit` returned by \link[rstan]{sampling}
}
\description{
Create an instance of the hierarchical SEIR Stan model incorporating
various data elements and sample model.
}
\examples{
stan_mod <- rstan::stan_model(system.file("stan",
  "hierarchical_SEIR_incidence_model.stan", package = "cr0eso"))
tmax <- 5
pop_size <- 100
dim(pop_size) <- c(1)
example_incidence <- matrix(c(1,1,2,3,2),ncol=1)
fit <- seir_model_fit(stan_model = stan_mod, tmax,1,example_incidence,pop_size)
}
\author{
Mike Irvine
}
