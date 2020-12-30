#' STAN helper functions including posterior sampling
#' SIR simulations and checks for stan model

# Load models
library(deSolve)
library(tidyverse)
library(rstan)
library(LaplacesDemon)
library(here)
library(conflicted)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# functions preferred
conflict_prefer("lag", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("extract", "rstan")
conflict_prefer("Blocks", "permute")
conflict_prefer("group_rows", "dplyr")
conflict_prefer("layout", "graphics")
conflict_prefer("lift", "purrr")
conflict_prefer("partial", "purrr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("chol2inv", "Matrix")
conflict_prefer("expand", "tidyr")
conflict_prefer("lmer", "lme4")

# Assign transmission and pathogen-induced death rates:
r0 = 2.5
sigma = 1/7
gamma = 1/7

# We will use the package deSolve to integrate, which requires certain data structures.
# Store parameters and initial values
# Parameters must be stored in a named list.
# dictionary: r0 - obvious, gamma - recovery rate, tau - time of start of interventions,
# zeta - strength of intervention
params <- list(r0 = r0,
               sigma = sigma,
               gamma = gamma,
               tau = 0,
               zeta = 0)





#' SEIR model for system of ODEs
#' See the 'ode' function documentation for further insights.
SEIR <- function(t, y, params) {
  with(as.list(c(params, y)), {
    S <- y[1]
    E <- y[2]
    I <- y[3]
    R <- y[4]

    dS = - r0 * gamma * S * I

    dE = r0 * gamma * S * I - sigma * E

    dI = sigma * E - gamma * I

    dR = gamma * I

    res <- c(dS,dE,dI,dR)
    list(res)
  })
}

#' Sample
#' @param out output from run_sim function
#' @param sample_days number of days sampled throughout the epidemic
#' @param sample_n number of host individuals sampled per day
#' @param bias Bias in sampling scheme bias = 0 means scheme purely random
#'             bias >0 means scheme postively biased (test seeking behaviour)
sampling_model <- function(out,tmax,sample_days = 20,
                           sample_n = 25,bias=0){


  # Choose which days the samples were taken.
  # Ideally this would be daily, but we all know that is difficult.
  sample_time = sort(sample(1:tmax, sample_days, replace=F))

  # Extract the "true" fraction of the population that is infected on each of the sampled days:
  sample_propinf = out[out$time %in% sample_time, 3]

  sample_propinf = invlogit(logit(sample_propinf) + bias)

  # Generate binomially distributed data.
  # So, on each day we sample a given number of people (sample_n), and measure how many are infected.
  # We expect binomially distributed error in this estimate, hence the random number generation.
  sample_y = rbinom(sample_days, sample_n, sample_propinf)
  return(list(y=sample_y,propinf=sample_propinf,
              time=sample_time))
}

#' Fr the SEIR function run similation and save as a data.frame
#' @param params parameter list
#' @param tmax maximum number of time-steps to simulate towards
run_sim <- function(params,tmax){
  E0 = 0.02    # initial fraction infected
  S0 = 1 - E0 # initial fraction susceptible
  I0 = 0
  R0 = 0

  # Initial conditions are stored in a vector
  inits <- c(S0, E0, I0, R0)

  # Create a time series over which to integrate.
  # Here we have an epidemic that is observed over tmax number of days (or weeks or etc).
  tmin = 0

  times = tmin:tmax
  # Run the integration:
  out <- ode(inits, times, SEIR, params, method="ode45")

  # Store the output in a data frame:
  out <- data.frame(out)
  colnames(out) <- c("time", "S", "E", "I", "R")

  return(out)
}

#' Fr the SEIR function run simulation with incidence and save as a data.frame
#' Includes decay of infectivity term to denote increasing intervention after
#' first case is detected
#' parameters:
#' * r0 - initial R_0 value
#' * gamma - recovery rate
#' * sigma rate of transition to infectious class
#' * tau - time of introduction of intervention
#' * zeta - strength of intervention
#' @param params parameter list
#' @param tmax maximum number of time-steps to simulate towards
#' @param E0 initial fraction infected
run_outbreak_sim <- function(params,tmax,E0 = 0.02){

  S0 = 1 - E0 # initial fraction susceptible
  I0 = 0
  R0 = 0
  Inc0 = 0

  # Initial conditions are stored in a vector
  inits <- c(S0, E0, I0, R0, Inc0)

  # define RHS
  RHS <- function(t, y, params) {
    with(as.list(c(params, y, t)), {
      S <- y[1]
      E <- y[2]
      I <- y[3]
      R <- y[4]
      Inc <- y[5]

      r <- r0*(t<=tau) +  (t>tau)*r0*exp(zeta*(tau-t))

      dS = - r * gamma * S * I

      dE = r * gamma * S * I - sigma * E

      dI = sigma * E - gamma * I

      dR = gamma * I

      dInc = sigma * E

      res <- c(dS, dE, dI, dR, dInc)
      list(res)
    })
  }

  # Create a time series over which to integrate.
  # Here we have an epidemic that is observed over tmax number of days (or weeks or etc).
  tmin = 0

  times = tmin:tmax
  # Run the integration:
  out <- ode(inits, times, RHS, params, method="ode45")

  # Store the output in a data frame:
  out <- data.frame(out)
  colnames(out) <- c("time", "S", "E", "I", "R", "Inc")

  out <- as_tibble(out) %>% mutate(I_diff = Inc - lag(Inc))

  return(out)
}

