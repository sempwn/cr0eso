



#' Extract model incidence and R0s from posterior
#' @author Mike Irvine
#' @note Naming convention throughout is snake case with prefix "hom_" to denote Hierarcical Outbreak Model
#' @param posts Object after calling extract of stan model object of hierarchical model
#' @importFrom tibble tibble as_tibble rowid_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate group_by arrange ungroup summarise
#' @export
hom_extract_posterior_draws <- function(posts){

  # extract posterior samples
  post_param <- tibble(r0=posts$predictive_r0,
                       sigma=posts$params[,2],
                       gamma=posts$params[,3],
                       n = posts$params[,4],
                       tau = posts$params[,5],
                       zeta = posts$zeta)

  # extract R0
  r0s <- as_tibble(posts$r0k) %>%
    pivot_longer(everything(),names_to="location",values_to="r0") %>%
    mutate(location = str_remove_all(location,"V"))

  # extract zeta (only applicable where zeta is multilevel)
  zetas <- NULL

  if("zetak" %in% names(posts)){
    zetas <- as_tibble(posts$zetak) %>%
      pivot_longer(everything(),names_to="location",values_to="zeta") %>%
      mutate(location = str_remove_all(location,"V"))
  }


  # extract cumulative incidence (cincidence) and then take diff to get incidence
  model_incidence <- as_tibble(posts$incidence) %>%
    rowid_to_column("run") %>%
    pivot_longer(-run,names_to="name",
                 values_to="cincidence") %>%
    separate(name,c("time","location")) %>%
    mutate(time = as.double(time)) %>%
    arrange(location,run,time) %>%
    group_by(location,run) %>%
    mutate(incidence = cincidence - lag(cincidence)) %>%
    ungroup()

  # simulate Poisson random variate and get credible intervals
  model_incidence <- model_incidence %>%
    mutate(rincidence = rpois(n(),incidence)) %>%
    group_by(time,location) %>%
    summarise(m = mean(rincidence,na.rm=TRUE),
              liqr = quantile(rincidence,0.25,na.rm=TRUE),
              uiqr = quantile(rincidence,0.75,na.rm=TRUE),
              lc = quantile(rincidence,0.025,na.rm=TRUE),
              uc = quantile(rincidence,0.975,na.rm=TRUE))

  return(list(r0=r0s,incidence=model_incidence,
              zeta=zetas,params=post_param))
}

#' Plot r0 by location as extracted from hierarchical outbreak model.
#' If extracted_posts is NULL then posts object is used to first extract the R0
#' and incidence draws from posterior (note this will take longer)
#' @param extracted_posts object returned by hom_extract_posterior_draws
#' @param posts Object after calling extract of stan model object of hierarchical model
#' @return list containing:
#'  * plot - ggplot object
#'  * table - tibble object of results
#' @examples
#'  mod # RSTAN fit object
#'  posts <- extract(mod) # extract posterior from model object
#'  extracted_posts <- hom_extract_posterior_draws(posts) # get object of incidence and r0
#'  result <- hom_plot_r0_by_location(extracted_posts=extracted_posts)
#'  # plot results
#'  show(result$plot)
#' @export
hom_plot_r0_by_location <- function(extracted_posts=NULL,posts=NULL){

  if(is.null(extracted_posts)){
    extracted_posts <- hom_extract_posterior_draws(posts)
  }


  r0s <- extracted_posts$r0
  params <- extracted_posts$params

  nlocs <- length(unique(r0s$location))
  ob_codes <- 1:nlocs

  # order locations by mean r0 (used to relevel location factor)
  loc_ordering <- r0s %>%
    group_by(location) %>%
    summarise(mean=mean(r0)) %>%
    arrange(mean) %>%
    pull(location)

  plot_data <- r0s %>%
    mutate(type="individual") %>%
    bind_rows(params %>%
                select(r0) %>%
                mutate(type="group") %>%
                mutate(location="Total")) %>%
    mutate(location = fct_rev(fct_relevel(location,
                                          c(loc_ordering,"Total")
    )))

  p <- plot_data %>% ggplot(aes(y=r0,x=location)) +
    geom_violin(alpha=0.5,aes(fill=type)) +

    stat_summary(fun.data = function(z){
      lc <- quantile(z,0.05)
      uc <- quantile(z,0.95)
      return(list(ymin=lc,ymax=uc))
    }, geom = "errorbar",size=1,width=0,
    color="gray40") +
    geom_point(size = 3, stat = 'summary',
               shape = 21, stroke = 2, fill="white",
               fun.y = function(x) mean(x, na.rm = TRUE)) +
    geom_hline(yintercept=1.0,linetype="dashed") +
    coord_flip() +
    theme_classic() +
    labs(y=expression(R[0])) +
    theme(legend.position = "none")

  r0_table <- plot_data %>%
    group_by(location) %>%
    summarise(m = quantile(r0,0.5),
              lc= quantile(r0,0.05),
              uc = quantile(r0,0.95)) %>%
    mutate_if(is.numeric,~round(.,2)) %>%
    mutate(r0 = glue::glue("{m} ({lc} - {uc})")) %>%
    mutate(ob_code = c(rev(ob_codes),"Total")) %>%
    select(ob_code,location,r0)
  list(plot=p,table=r0_table)
}

#' Plot zeta by location as extracted from hierarchical outbreak model.
#' If extracted_posts is NULL then posts object is used to first extract the zeta
#' and incidence draws from posterior (note this will take longer)
#' @param extracted_posts object returned by hom_extract_posterior_draws
#' @param posts Object after calling extract of stan model object of hierarchical model
#' @return list containing:
#'  * plot - ggplot object
#'  * table - tibble object of results
#' @examples
#'  mod # RSTAN fit object
#'  posts <- extract(mod) # extract posterior from model object
#'  extracted_posts <- hom_extract_posterior_draws(posts) # get object of incidence and zeta
#'  result <- hom_plot_zeta_by_location(extracted_posts=extracted_posts)
#'  # plot results
#'  show(result$plot)
#' @export
hom_plot_zeta_by_location <- function(extracted_posts=NULL,posts=NULL){

  if(is.null(extracted_posts)){
    extracted_posts <- hom_extract_posterior_draws(posts)
  }

  zetas <- extracted_posts$zeta
  params <- extracted_posts$params

  nlocs <- length(unique(zetas$location))
  ob_codes <- 1:nlocs

  # order locations by mean zeta (used to relevel location factor)
  loc_ordering <- zetas %>%
    group_by(location) %>%
    summarise(mean=mean(zeta)) %>%
    arrange(mean) %>%
    pull(location)

  plot_data <- zetas %>%
    mutate(type="individual") %>%
    bind_rows(params %>%
                select(zeta) %>%
                mutate(type="group") %>%
                mutate(location="Total")) %>%
    mutate(location = fct_rev(fct_relevel(location,
                                          c(loc_ordering,"Total")
    )))

  p <- plot_data %>% ggplot(aes(y=zeta,x=location)) +
    geom_violin(alpha=0.5,aes(fill=type)) +

    stat_summary(fun.data = function(z){
      lc <- quantile(z,0.05)
      uc <- quantile(z,0.95)
      return(list(ymin=lc,ymax=uc))
    }, geom = "errorbar",size=1,width=0,
    color="gray40") +
    geom_point(size = 3, stat = 'summary',
               shape = 21, stroke = 2, fill="white",
               fun.y = function(x) mean(x, na.rm = TRUE)) +
    coord_flip() +
    theme_classic() +
    labs(y=expression(paste("Intervention strength (",zeta,")"))) +
    theme(legend.position = "none")

  zeta_table <- plot_data %>%
    group_by(location) %>%
    summarise(m = quantile(zeta,0.5),
              lc= quantile(zeta,0.05),
              uc = quantile(zeta,0.95)) %>%
    mutate_if(is.numeric,~round(.,2)) %>%
    mutate(zeta = glue::glue("{m} ({lc} - {uc})")) %>%
    mutate(ob_code = c(rev(ob_codes),"Total")) %>%
    select(ob_code,location,zeta)
  list(plot=p,table=zeta_table)
}


#' Plot incidence by location as extracted from hierarchical outbreak model.
#' If extracted_posts is NULL then posts object is used to first extract the
#' incidence draws from posterior (note this will take longer)
#' @param extracted_posts object returned by hom_extract_posterior_draws
#' @param posts Object after calling extract of stan model object of hierarchical model
#' @param outbreak_cases matrix of outbreak case data
#' @param end_time time in days to plot up until
#' @return list containing:
#'  * plot - ggplot object
#'  * table - tibble object of results
#' @examples
#'  mod # RSTAN fit object
#'  outbreak_cases # data of outbreaks in matrix form
#'  posts <- extract(mod) # extract posterior from model object
#'  extracted_posts <- hom_extract_posterior_draws(posts) # get object of incidence and r0
#'  result <- hom_plot_incidence_by_location(extracted_posts=extracted_posts,
#'                                           outbreak_cases=outbreak_cases)
#'  # plot results
#'  show(result$plot)
#' @export
hom_plot_incidence_by_location <- function(extracted_posts=NULL,posts=NULL,
                                    outbreak_cases=NULL, end_time=60){

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
    group_by(location) %>%
    summarise(mean=mean(r0)) %>%
    arrange(mean) %>%
    mutate(label = glue::glue("Location: {location}, R0: {round(mean,2)}"))

  # get ordering of locations by R0
  loc_ordering <- loc_labels %>%
    pull(location)

  # transform outbreak case data for plotting
  cases_data <- as_tibble(outbreak_cases)
  names(cases_data) <- as.character(1:nlocs)
  cases_data <- cases_data %>%
    rowid_to_column("time") %>%
    pivot_longer(-time,names_to="location",values_to="data_incidence")

  # join posterior and data for plotting
  plot_data <- incidence %>%
    inner_join(cases_data,by=c("time","location")) %>%
    mutate(time = as.double(time)) %>%
    mutate(location = fct_relevel(location,loc_ordering)) %>%
    mutate(
      location = fct_relabel(location,
                             function(x){
                               loc_labels %>%
                                 filter(location==x) %>%
                                 pull(label)
                             })
    ) %>%
    filter(time < end_time)

  p <- plot_data %>%
    ggplot(aes(x=time,group=location)) +
    geom_ribbon(aes(ymin=lc,ymax=uc),fill="steelblue",alpha=0.3) +
    geom_ribbon(aes(ymin=liqr,ymax=uiqr),fill="steelblue",alpha=0.3) +
    geom_line(aes(y=m)) +
    geom_point(aes(y=data_incidence),fill="white",size=1.5) +
    facet_wrap(vars(location),ncol=3,scales="free_y") +
    theme_classic() +
    theme(strip.text = element_text(size=11)) +
    labs(y="Cases",x="Time since initial case (days)")

  list(plot=p,data=plot_data)
}
