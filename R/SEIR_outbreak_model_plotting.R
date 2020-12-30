#' Outbreak model plots
#' @author Mike Irvine
#' 

library(tidyverse)
library(here)

#' Plot posterior and prior density comparison
#' @param prior_param dataframe of prior samples where each row is a sample and each column is a parameter
#' @param post_param dataframe of posterior samples in same format as prior_param
#' @param true_vals list of true values if exist. List should contain parameters and their corresponding values
prior_post_plot <- function(prior_param,post_param,true_vals = NULL){
  plot_data <- prior_param %>%
    mutate(type="prior") %>%
    slice(1000:n()) %>% # remove burn-in
    bind_rows(
      post_param %>%
        mutate(type="posterior") %>%
        slice(1000:n())      
    ) %>%
    pivot_longer(-type,names_to="param",values_to = "value") 
  
  p <- plot_data %>% 
    ggplot() + 
    geom_density(aes(x=value,fill=type,color=type),alpha=0.5) 
  
  if(!is.null(true_vals)){
    val_data <- true_vals %>% 
      as_tibble() %>% 
      pivot_longer(everything(),
                   names_to="param",
                   values_to = "true_value")
    
    p <- p +
      geom_vline(data=val_data, aes(xintercept=true_value)) +
      facet_wrap(vars(param), scales="free")
  }
  
  p <- p +
    facet_wrap(vars(param), scales="free")
  
  
}

#' Plot model predictions across the sampling time period with data
#' @param post_data posterior from STAN output in samples x time format
#' @param df_sample dataframe of data must contain time and y as data output
plot_from_sample_post <- function(post_data,df_sample){
  mod_median = apply(post_data, 2, median, na.rm = TRUE)
  mod_low = apply(post_data, 2, quantile, probs=c(0.025), na.rm = TRUE)
  mod_high = apply(post_data, 2, quantile, probs=c(0.975), na.rm = TRUE)
  mod_iqr_low = apply(post_data, 2, quantile, probs=c(0.25), na.rm = TRUE)
  mod_iqr_high = apply(post_data, 2, quantile, probs=c(0.75), na.rm = TRUE)
  mod_time = 1:length(mod_median)
  
  df_fit = tibble(mod_median, mod_low, mod_high,
                  mod_iqr_low, mod_iqr_high,
                  mod_time)
  
  # Plot the synthetic data with the model predictions
  # Median and 95% Credible Interval
  
  ggplot() +
    
    # Error in integration:
    
    geom_ribbon(data = df_fit, aes(x=lubridate::ymd("2020-03-01") + lubridate::days(mod_time),
                                   ymin=mod_low,
                                   ymax=mod_high, 
                                   fill = "PPC"),alpha=0.3) +
    geom_ribbon(data = df_fit, aes(x=lubridate::ymd("2020-03-01") + lubridate::days(mod_time),
                                   ymin=mod_iqr_low,
                                   ymax=mod_iqr_high, 
                                   fill = "PPC"),alpha=0.3) +
    geom_line(data = df_fit, aes(x=lubridate::ymd("2020-03-01") + lubridate::days(mod_time),
                                 y=mod_median, fill = "PPC")) + 
    geom_point(data=df_sample, aes(x=lubridate::ymd("2020-03-01") + lubridate::days(time),
                                   y=y),
               col="black", shape = 19, size = 1.5) +
    # Aesthetics
    labs(x = "Time (days)", y = "Proportion Infected") + 
    theme_classic() + 
    theme(axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"))
  
}

#' Plot model predictions across the sampling time period with data
#' @param posts posterior from STAN output
#' @param samples data samples of number infected
#' @param index the index for output from the STAN ODE model
#' @param incidence boolean whether plotting incidence or plotting a sample
plot_posterior_data <- function(posts,samples,index=3,
                                incidence = FALSE){
  
  # Combine into two data frames for plotting
  if(incidence){
    df_sample = as_tibble(samples) %>%
      mutate(y = data_y)
  } else{
    df_sample = as_tibble(samples) %>%
      mutate(y = y/sample_n)
  }
  plot_from_sample_post(posts$fake_I[,,index],df_sample)
}
