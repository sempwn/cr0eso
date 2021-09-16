



#' Extract model incidence and R0s from posterior
#' @note Naming convention throughout is snake case with prefix "hom_" to denote Hierarcical Outbreak Model
#' @param posts Object after calling extract of stan model object of hierarchical model
#' @param location_labels tibble with columns `site` and `location`. Will be used to
#'   relabel numeric location to values in `site` column
#' @importFrom stats quantile rpois
#' @importFrom magrittr %>%
#' @export
hom_extract_posterior_draws <- function(posts,
                                        location_labels = NULL ){

  location <- run <- name <- time <- cincidence <- incidence <- rincidence <- NULL


  # extract posterior samples
  post_param <- tibble::tibble(r0=posts$predictive_r0,
                       sigma=posts$params[,2],
                       gamma=posts$params[,3],
                       n = posts$params[,4],
                       tau = posts$params[,5],
                       zeta = posts$zeta)

  # extract R0
  r0s <- tibble::as_tibble(posts$r0k) %>%
    tidyr::pivot_longer(dplyr::everything(),names_to="location",values_to="r0") %>%
    dplyr::mutate(location = stringr::str_remove_all(location,"V"))

  # extract zeta (only applicable where zeta is multilevel)
  zetas <- NULL

  if("zetak" %in% names(posts)){
    zetas <- tibble::as_tibble(posts$zetak) %>%
      tidyr::pivot_longer(dplyr::everything(),names_to="location",values_to="zeta") %>%
      dplyr::mutate(location = stringr::str_remove_all(location,"V"))
  }


  # extract cumulative incidence (cincidence) and then take diff to get incidence
  model_incidence <- tibble::as_tibble(posts$incidence) %>%
    tibble::rowid_to_column("run") %>%
    tidyr::pivot_longer(-run,names_to="name",
                 values_to="cincidence") %>%
    tidyr::separate(name,c("time","location")) %>%
    dplyr::mutate(time = as.double(time)) %>%
    dplyr::arrange(location,run,time) %>%
    dplyr::group_by(location,run) %>%
    dplyr::mutate(incidence = cincidence - dplyr::lag(cincidence)) %>%
    dplyr::ungroup()

  # simulate Poisson random variate and get credible intervals
  model_incidence <- model_incidence %>%
    dplyr::mutate(rincidence = rpois(dplyr::n(),incidence)) %>%
    dplyr::group_by(time,location) %>%
    dplyr::summarise(m = mean(rincidence,na.rm=TRUE),
              liqr = quantile(rincidence,0.25,na.rm=TRUE),
              uiqr = quantile(rincidence,0.75,na.rm=TRUE),
              lc = quantile(rincidence,0.025,na.rm=TRUE),
              uc = quantile(rincidence,0.975,na.rm=TRUE))

  # update labels if provided
  if(!is.null(location_labels)){
    r0s <- r0s %>%
      replace_location_values(location_labels)

    if("zetak" %in% names(posts)){
      zetas <- zetas %>%
        replace_location_values(location_labels)
    }

    model_incidence <- model_incidence %>%
      replace_location_values(location_labels)

  }

  return(list(r0=r0s,incidence=model_incidence,
              zeta=zetas,params=post_param))
}

#' Plot r0 by location as extracted from hierarchical outbreak model.
#' @description If extracted_posts is NULL then posts object is used to first extract the R0
#' and incidence draws from posterior (note this will take longer)
#' @note Naming convention throughout is snake case with prefix "hom_" to denote Hierarchical Outbreak Model
#' @param extracted_posts object returned by hom_extract_posterior_draws
#' @param posts Object after calling extract of stan model object of hierarchical model
#' @param sort TRUE/FALSE. If TRUE (default), locations will be ordered by mean R0
#' @importFrom stats quantile
#' @return list containing:
#'  * plot - ggplot object
#'  * table - tibble object of results
#' @examples
#' \dontrun{
#'  mod # RSTAN fit object
#'  posts <- rstan::extract(mod) # extract posterior from model object
#'  extracted_posts <- hom_extract_posterior_draws(posts) # get object of incidence and r0
#'  result <- hom_plot_r0_by_location(extracted_posts=extracted_posts)
#'  # plot results
#'  show(result$plot)
#'  }
#' @export
hom_plot_r0_by_location <- function(extracted_posts=NULL,posts=NULL,sort=TRUE){

  r0 <- type <- ob_code <- location <- NULL

  if(is.null(extracted_posts)){
    extracted_posts <- hom_extract_posterior_draws(posts)
  }


  r0s <- extracted_posts$r0
  params <- extracted_posts$params

  nlocs <- length(unique(r0s$location))
  ob_codes <- 1:nlocs

  # order locations by mean r0 (used to relevel location factor)
  loc_ordering <- r0s %>%
    dplyr::group_by(location) %>%
    dplyr::summarise(mean=mean(r0)) %>%
    dplyr::arrange(if(sort==TRUE){mean}else{location}) %>%
    dplyr::pull(location)

  plot_data <- r0s %>%
    dplyr::mutate(type="individual") %>%
    dplyr::bind_rows(params %>%
                dplyr::select(r0) %>%
                dplyr::mutate(type="group") %>%
                dplyr::mutate(location="Total")) %>%
    dplyr::mutate(location = forcats::fct_rev(forcats::fct_relevel(location,
                                          c(loc_ordering,"Total")
    )))

  p <- plot_data %>%
    ggplot2::ggplot(ggplot2::aes(y=r0,x=location)) +
    ggplot2::geom_violin(alpha=0.5,ggplot2::aes(fill=type)) +

    ggplot2::stat_summary(fun.data = function(z){
      lc <- quantile(z,0.05)
      uc <- quantile(z,0.95)
      return(list(ymin=lc,ymax=uc))
    }, geom = "errorbar",size=1,width=0,
    color="gray40") +
    ggplot2::geom_point(size = 3, stat = 'summary',
               shape = 21, stroke = 2, fill="white",
               fun.y = function(x) mean(x, na.rm = TRUE)) +
    ggplot2::geom_hline(yintercept=1.0,linetype="dashed") +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::labs(y=expression(R[0])) +
    ggplot2::theme(legend.position = "none")

  r0_table <- plot_data %>%
    dplyr::group_by(location) %>%
    dplyr::summarise(m = quantile(r0,0.5),
              lc= quantile(r0,0.05),
              uc = quantile(r0,0.95)) %>%
    dplyr::mutate_if(is.numeric,~round(.,2)) %>%
    dplyr::mutate(r0 = glue::glue("{m} ({lc} - {uc})")) %>%
    #dplyr::mutate(ob_code = c(rev(ob_codes),"Total")) %>%
    dplyr::select(location,r0)
  list(plot=p,table=r0_table)
}

#' Plot zeta by location as extracted from hierarchical outbreak model.
#' @description If extracted_posts is NULL then posts object is used to first extract the zeta
#' and incidence draws from posterior (note this will take longer)
#' @note Naming convention throughout is snake case with prefix "hom_" to denote Hierarcical Outbreak Model
#' @param extracted_posts object returned by hom_extract_posterior_draws
#' @param posts Object after calling extract of stan model object of hierarchical model
#' @param sort TRUE/FALSE. If TRUE (default), locations will be ordered by mean R0
#' @importFrom stats quantile
#' @importFrom magrittr %>%
#' @return list containing:
#'  * plot - ggplot object
#'  * table - tibble object of results
#' @examples
#' \dontrun{
#'  mod # RSTAN fit object
#'  posts <- rstan::extract(mod) # extract posterior from model object
#'  extracted_posts <- hom_extract_posterior_draws(posts) # get object of incidence and zeta
#'  result <- hom_plot_zeta_by_location(extracted_posts=extracted_posts)
#'  # plot results
#'  show(result$plot)
#'  }
#' @export
hom_plot_zeta_by_location <- function(extracted_posts=NULL,posts=NULL,sort=TRUE){

  location <- zeta <- type <- ob_code <- NULL

  if(is.null(extracted_posts)){
    extracted_posts <- hom_extract_posterior_draws(posts)
  }

  zetas <- extracted_posts$zeta
  params <- extracted_posts$params

  nlocs <- length(unique(zetas$location))
  ob_codes <- 1:nlocs

  # order locations by mean zeta (used to relevel location factor)
  loc_ordering <- zetas %>%
    dplyr::group_by(location) %>%
    dplyr::summarise(mean=mean(zeta)) %>%
    dplyr::arrange(if(sort==TRUE){mean}else{location}) %>%
    dplyr::pull(location)

  plot_data <- zetas %>%
    dplyr::mutate(type="individual") %>%
    dplyr::bind_rows(zetas %>%
                       dplyr::select(zeta) %>%
                       dplyr::mutate(type="group") %>%
                       dplyr::mutate(location="Total")) %>%
    dplyr::mutate(location = forcats::fct_rev(forcats::fct_relevel(location,
                                          c(loc_ordering,"Total")
    )))

  p <- plot_data %>%
    ggplot2::ggplot(ggplot2::aes(y=zeta,x=location)) +
    ggplot2::geom_violin(alpha=0.5,ggplot2::aes(fill=type)) +

    ggplot2::stat_summary(fun.data = function(z){
      lc <- quantile(z,0.05)
      uc <- quantile(z,0.95)
      return(list(ymin=lc,ymax=uc))
    }, geom = "errorbar",size=1,width=0,
    color="gray40") +
    ggplot2::geom_point(size = 3, stat = 'summary',
               shape = 21, stroke = 2, fill="white",
               fun.y = function(x) mean(x, na.rm = TRUE)) +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::labs(y=expression(paste("Intervention strength (",zeta,")"))) +
    ggplot2::theme(legend.position = "none")

  zeta_table <- plot_data %>%
    dplyr::group_by(location) %>%
    dplyr::summarise(m = quantile(zeta,0.5),
              lc= quantile(zeta,0.05),
              uc = quantile(zeta,0.95)) %>%
    dplyr::mutate_if(is.numeric,~round(.,2)) %>%
    dplyr::mutate(zeta = glue::glue("{m} ({lc} - {uc})")) %>%
    #dplyr::mutate(ob_code = c(rev(ob_codes),"Total")) %>%
    dplyr::select(ob_code,location,zeta)
  list(plot=p,table=zeta_table)
}


#' Plot the critical time for each location
#' @description Critical time defined as time in days until the effective
#'  R number reduces to one or below. Plot also gives 95 credible interval and
#'  attack rate for each location.
#' @note Naming convention throughout is snake case with prefix "hom_" to denote Hierarcical Outbreak Model
#' @param posts Object after calling extract of stan model object of hierarchical model
#' @param locations vector of location names
#' @param num_cases vector of number of cases for each location
#' @param capacities vector of sizes for each location
#' @param location_labels tibble with columns `site` and `location`. Will be used to
#'   relabel numeric location to values in `site` column
#' @importFrom stats quantile rpois
#' @importFrom magrittr %>%
#' @export
hom_plot_critical_times_by_location <- function(posts,
                                                locations,
                                                num_cases,
                                                capacities,
                                                location_labels = NULL){


  attack_rate <- location <- critical_time <- NULL



  # Calculate attack rates
  ar_table <- data.frame("location" = as.character(locations),
                         "attack_rate" = 100*num_cases/capacities)

  # get critical times and combine
  plot_data <- create_critical_times(posts) %>%
    dplyr::group_by(location) %>%
    dplyr::summarise(m = stats::median(critical_time),
                     lc = stats::quantile(critical_time,0.05),
                     uc = stats::quantile(critical_time,0.95)) %>%
    dplyr::left_join(ar_table,by="location")

  if(!is.null(location_labels)){
    plot_data <- plot_data %>%
      replace_location_values(location_labels)
  }



  # plot data
  plot_data %>% ggplot2::ggplot(ggplot2::aes(x=location, y=m)) +
    ggplot2::geom_point(aes(colour = attack_rate/100), size=2.5, stat="identity") +
    ggplot2::geom_errorbar(aes(ymin = lc, ymax = uc, colour = attack_rate/100), width = 0.6, size=0.6) +
    ggplot2::scale_colour_continuous(name = "Attack rate",  labels = scales::percent,
                                     type="gradient", low="#DF536B",
                                     high="black") +
    ggplot2::ylab("Critical time (days)") +
    ggplot2::xlab("Location") +
    ggplot2::theme_classic()
}


#' Plot incidence by location as extracted from hierarchical outbreak model.
#' @description If extracted_posts is NULL then posts object is used to first extract the
#' incidence draws from posterior (note this will take longer)
#' @note Naming convention throughout is snake case with prefix "hom_" to denote Hierarcical Outbreak Model
#' @param extracted_posts object returned by hom_extract_posterior_draws
#' @param posts Object after calling extract of stan model object of hierarchical model
#' @param outbreak_cases matrix of outbreak case data
#' @param end_time time in days to plot up until
#' @param location_labels tibble with columns `site` and `location`. Will be used to
#'   relabel numeric location to values in `site` column
#' @param sort TRUE/FALSE. If TRUE (default), locations will be ordered by mean R0
#' @importFrom magrittr %>%
#' @return list containing:
#'  * plot - ggplot object
#'  * table - tibble object of results
#' @examples
#' \dontrun{
#'  mod # RSTAN fit object
#'  outbreak_cases <- matrix(c(0,0,1,2,3,0,0,2,4,5),ncol=2) # data of outbreaks in matrix form
#'  posts <- rstan::extract(mod) # extract posterior from model object
#'  extracted_posts <- hom_extract_posterior_draws(posts) # get object of incidence and r0
#'  result <- hom_plot_incidence_by_location(extracted_posts=extracted_posts,
#'                                           outbreak_cases=outbreak_cases)
#'  # plot results
#'  show(result$plot)
#'  }
#' @export
hom_plot_incidence_by_location <- function(extracted_posts=NULL,posts=NULL,
                                    outbreak_cases=NULL, end_time=60,
                                    location_labels = NULL,sort=TRUE){

  r0 <- time <- location <- label <- lc <- uc <- liqr <- uiqr <- m <- data_incidence <- NULL

  if(is.null(extracted_posts)){
    extracted_posts <- hom_extract_posterior_draws(posts)
  }

  incidence <- extracted_posts$incidence
  r0s <- extracted_posts$r0
  params <- extracted_posts$params
  nlocs <- length(unique(r0s$location))

  # create factor reference tibble to
  # relabel factors into their mean R0
  loc_labels <- r0s %>%
    dplyr::group_by(location) %>%
    dplyr::summarise(mean=mean(r0)) %>%
    dplyr::arrange(if(sort==TRUE){mean}else{location}) %>%
    dplyr::mutate(label = glue::glue("Location: {location}, R0: {round(mean,2)}"))

  # get ordering of locations by R0
  loc_ordering <- loc_labels %>%
    dplyr::pull(location)

  # transform outbreak case data for plotting
  cases_data <- tibble::as_tibble(outbreak_cases)
  names(cases_data) <- as.character(1:nlocs)
  cases_data <- cases_data %>%
    tibble::rowid_to_column("time") %>%
    tidyr::pivot_longer(-time,names_to="location",values_to="data_incidence")

  if(!is.null(location_labels)){
    cases_data <- cases_data %>%
      replace_location_values(location_labels)
  }

  # join posterior and data for plotting
  plot_data <- incidence %>%
    dplyr::inner_join(cases_data,by=c("time","location")) %>%
    dplyr::mutate(time = as.double(time)) %>%
    dplyr::mutate(location = forcats::fct_relevel(location,loc_ordering)) %>%
    dplyr::mutate(
      location = forcats::fct_relabel(location,
                             function(x){
                               loc_labels %>%
                                 dplyr::filter(location==x) %>%
                                 dplyr::pull(label)
                             })
    ) %>%
    dplyr::filter(time < end_time)

  p <- plot_data %>%
    ggplot2::ggplot(ggplot2::aes(x=time, group=location)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lc,ymax=uc),fill="steelblue",alpha=0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=liqr,ymax=uiqr),fill="steelblue",alpha=0.3) +
    ggplot2::geom_line(ggplot2::aes(y=m)) +
    ggplot2::geom_point(ggplot2::aes(y=data_incidence),fill="white",size=1.5) +
    ggplot2::facet_wrap(ggplot2::vars(location),ncol=3,scales="free_y") +
    ggplot2::theme_classic() +
    ggplot2::theme(strip.text = ggplot2::element_text(size=11)) +
    ggplot2::labs(y="Cases",x="Time since initial case (days)")

  list(plot=p,data=plot_data)
}


#' plot incidence and counterfactual scenario
#' @description on same plot show counterfactual and incidence posterior predictive distributions
#' across time and faceted by location. Shaded regions represent the 50% and 90% CrI, with median of
#' draws shown as lines.
#' @note Naming convention throughout is snake case with prefix "hom_" to denote Hierarcical Outbreak Model
#' @param fit output of [seir_model_fit]
#' @param outbreak_cases matrix of outbreak case data
#' @param location_labels tibble with columns `site` and `location`. Will be used to
#'   relabel numeric location to values in `site` column
#' @importFrom magrittr %>%
#' @importFrom stats quantile median
#' @return list containing:
#'  * plot - ggplot object
#'  * table - tibble object of results
#'  * violin_plot - ggplot object
#' @examples
#' \dontrun{
#'  tmax <- 5
#'  pop_size <- 100
#'  dim(pop_size) <- c(1)
#'  example_incidence <- matrix(c(1,1,2,3,2),ncol=1)
#'  fit <- seir_model_fit(tmax,1,example_incidence,pop_size)
#'  result <- hom_plot_incidence_by_location(fit)
#'  # plot results
#'  show(result$plot)
#'  }
#' @export
hom_plot_counterfactual_by_location <- function(fit, outbreak_cases = NULL,
                                                location_labels = NULL){

  location <- counterfactual_cases <- scenario <- .draw <- cases <- baseline <- NULL
  averted <- proportion_averted <- m <- liqr <- uiqr <- lc <- uc <- data_incidence <- NULL
  time <- NULL

  # extract posterior predictive distribution
  pp_cases <- fit$model %>%
    tidybayes::spread_draws(pp_cases[time,location]) %>%
    dplyr::filter(time > 1)

  # extract posterior counterfactual distribution
  counterfactual <- fit$model %>%
    tidybayes::spread_draws(counterfactual_cases[time,location]) %>%
      dplyr::filter(time > 1)

  if(!is.null(location_labels)){
    counterfactual <- counterfactual %>%
      replace_location_values(location_labels)

    pp_cases <- pp_cases %>%
      replace_location_values(location_labels)

  }

  # transform outbreak case data for plotting
  if(!is.null(outbreak_cases)){
    cases_data <- tibble::as_tibble(outbreak_cases)
    names(cases_data) <- as.character(1:ncol(cases_data))
    cases_data <- cases_data %>%
      tibble::rowid_to_column("time") %>%
      tidyr::pivot_longer(-time,names_to="location",values_to="data_incidence")

    if(!is.null(location_labels)){
      cases_data <- cases_data %>%
        replace_location_values(location_labels)
    }
  }

  plot_data <- dplyr::inner_join(pp_cases,counterfactual,
                                 by=c("time","location",
                                      ".chain",".iteration",".draw")) %>%
    dplyr::mutate(location = as.character(location)) %>%
    tidyr::pivot_longer(c("pp_cases","counterfactual_cases"),
                        names_to="scenario",
                        values_to="cases") %>%
    dplyr::mutate(scenario = dplyr::recode(scenario,
                                           "pp_cases" = "baseline",
                                           "counterfactual_cases" = "counterfactual"))

  # create counterfactual summary table
  counterfactual_table <- plot_data %>%
    # note draw is unique code for each draw from posterior
    dplyr::group_by(location,scenario,.draw) %>%
    dplyr::summarise(total_cases = sum(cases)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from="scenario",values_from="total_cases")

  # create total row by first summing across all locations and then summarising
  total_counterfactual_table <- counterfactual_table %>%
    dplyr::group_by(.draw) %>%
    dplyr::summarise(counterfactual = sum(counterfactual),
                     baseline = sum(baseline)
                     ) %>%
    dplyr::mutate(averted = counterfactual - baseline,
                  proportion_averted = averted/counterfactual
    ) %>%
    dplyr::summarise(baseline = summary_m_95cri(baseline),
                     counterfactual = summary_m_95cri(counterfactual),
                     averted = summary_m_95cri(averted),
                     proportion_averted = summary_m_95cri(proportion_averted)
    ) %>%
    dplyr::mutate(location = "total")

  # get data to create the violin plots of proportion averted
  violin_plot_data <- counterfactual_table %>%
    dplyr::mutate(averted = counterfactual - baseline,
                  proportion_averted = averted/counterfactual
    )

  counterfactual_table <- counterfactual_table %>%
    dplyr::mutate(averted = counterfactual - baseline,
                  proportion_averted = averted/counterfactual
    ) %>%
    dplyr::group_by(location) %>%
    dplyr::summarise(baseline = summary_m_95cri(baseline),
                     counterfactual = summary_m_95cri(counterfactual),
                     averted = summary_m_95cri(averted),
                     proportion_averted = summary_m_95cri(proportion_averted)
                     ) %>%
    dplyr::mutate(location = as.character(location)) %>%
    dplyr::bind_rows(total_counterfactual_table)


  # create 5,25,50,75,95 percentile by time
  plot_data <- plot_data %>%
    dplyr::group_by(time,location,scenario) %>%
    dplyr::summarise(m = median(cases,na.rm=TRUE),
              liqr = quantile(cases,0.25,na.rm=TRUE),
              uiqr = quantile(cases,0.75,na.rm=TRUE),
              lc = quantile(cases,0.05,na.rm=TRUE),
              uc = quantile(cases,0.95,na.rm=TRUE))


  # if including data then add to plot_data object
  if(!is.null(outbreak_cases)){
    plot_data <- plot_data %>%
      dplyr::inner_join(cases_data,by=c("time","location"))
  }

  p <- plot_data %>%
    ggplot2::ggplot(ggplot2::aes(x=time,y=m,fill=scenario)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=liqr,ymax=uiqr),alpha=0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lc,ymax=uc),alpha=0.3) +
    ggplot2::geom_line(ggplot2::aes(color=scenario))

  # if including data then add to plot
  if(!is.null(outbreak_cases)){
    p <- p +
      ggplot2::geom_point(ggplot2::aes(y=data_incidence),
                          fill="white",size=1.5)
  }

  p <- p +
    ggplot2::facet_wrap(ggplot2::vars(location), ncol=3,scales="free_y") +
    ggplot2::theme_classic() +
    ggplot2::labs(y="Cases",x="Time since initial case (days)")

  violin_plot <- violin_plot_data %>%
    ggplot2::ggplot(ggplot2::aes(x=location,y=proportion_averted)) +
    ggplot2::geom_violin(fill = "#3366FF",
                         draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::coord_cartesian(ylim=c(0,1), expand=FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_classic() +
    ggplot2::labs(y="Proportion cases averted",
                  x="Location")

  list(plot=p,table=counterfactual_table,
       violin_plot = violin_plot)
}


#' convenience function replace location labels with site variables from
#' location variable tables
#' @param d dataframe
#' @param location_labels dataframe with columns `location` and `site`
#' @noRd
replace_location_values <- function(d,location_labels){
  location <- NULL
  d %>%
    dplyr::ungroup() %>%
    dplyr::mutate(location = as.character(location)) %>%
    dplyr::left_join(location_labels,by="location") %>%
    dplyr::select(-location) %>%
    dplyr::rename("location" = "site")
}


