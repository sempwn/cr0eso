
#' summarise a model output by location
#' @param x location based output tibble
#' @param post_param overall (Total) for output tibble
#' @param output str name of output to create column for
#' @importFrom rlang :=
#' @returns tibble with columns location and "output"
#' @noRd
summarise_by_location <- function(x,post_param,output="r0"){
  location <- NULL

  nlocs <- length(unique(x$location))

  plot_data <- x %>%
    dplyr::mutate(type="individual") %>%
    dplyr::bind_rows(post_param %>%
                       dplyr::select({{output}}) %>%
                       dplyr::mutate(type="group") %>%
                       dplyr::mutate(location="Total")) %>%
    dplyr::mutate(location = forcats::fct_rev(forcats::fct_relevel(location,
                                                                   c("Total",as.character(1:nlocs))
    )))

  # convert to summary table
  plot_data <- plot_data %>%
    dplyr::group_by(location) %>%
    dplyr::summarise(m = quantile(!!rlang::sym(output),0.5, na.rm=TRUE),
                     lc= quantile(!!rlang::sym(output),0.05, na.rm=TRUE),
                     uc = quantile(!!rlang::sym(output),0.95, na.rm=TRUE)) %>%
    dplyr::mutate_if(is.numeric,~round(.,2)) %>%
    dplyr::mutate(!!output := glue::glue("{m} ({lc} - {uc})")) %>%
    dplyr::select(location,{{output}})

  return(plot_data)
}

#' Create time to hit r0 summary
#' @noRd
create_critical_times <- function(posts){

  ID <- r0 <- zeta <- location <- NULL

  r0s <- tibble::as_tibble(posts$r0k) %>%
    tibble::rowid_to_column("ID") %>%
    tidyr::pivot_longer(-ID,names_to="location",values_to="r0")

  zetas <- tibble::as_tibble(posts$zetak) %>%
    tibble::rowid_to_column("ID") %>%
    tidyr::pivot_longer(-ID,names_to="location",values_to="zeta")

  # create time to hit R0 summary
  crit_times <-
    dplyr::inner_join(r0s,zetas,by=c("ID","location")) %>%
    # critical time is time to reach R0 in days
    dplyr::mutate(critical_time = dplyr::if_else(r0<=1,0,log(r0)/zeta)) %>%
    dplyr::mutate(location = stringr::str_remove_all(location,"V"))

  return(crit_times)
}

#' Create comparison table for multiple posteriors
#' @description Create a comparison table of posteriors generated from the same
#' data. For example you may wish to compare the R0 where an intervention is
#' assumed compared to another model where no intervention is assumed
#' @param ... named list where each item is output of [covid_fit_seir]
#' @returns tibble of model results
#' @importFrom rlang :=
#' @examples
#' \dontrun{
#' create_pub_tables(mod)
#' create_pub_tables(intervention = mod,no_intervention= mod2)
#' }
#' @export
create_pub_tables <- function(...){

  location <- NULL

  # create list of cr0eso outputs
  mod_list <- list(...)

  # create placeholder list to fill with results for each output
  res_list <- list()

  # create index
  i <- 1
  for(model_name in names(mod_list)){

    # get model posterior
    mod <- mod_list[[model_name]]
    posts <- rstan::extract(mod$model)

    # get posterior parameters
    mod_post <- tibble::tibble(r0=posts$r0,
                               zeta = posts$zeta)

    # get r0s by location
    mod_r0s <- tibble::as_tibble(posts$r0k) %>%
      tidyr::pivot_longer(tidyselect::everything(),names_to="location",values_to="r0") %>%
      dplyr::mutate(location = stringr::str_remove_all(location,"V"))

    mod_zetas <- tibble::as_tibble(posts$zetak) %>%
      tidyr::pivot_longer(tidyselect::everything(),names_to="location",values_to="zeta") %>%
      dplyr::mutate(location = stringr::str_remove_all(location,"V"))

    # time to r=1
    mod_times <- create_critical_times(posts)

    # create labels for columns
    r0_label <- paste("r0",model_name)
    zeta_label <- paste("zeta",model_name)
    critical_time_label <- paste("critical_time",model_name)


    res <- summarise_by_location(mod_r0s, mod_post) %>%
      dplyr::rename(!!r0_label := "r0") %>%
      dplyr::inner_join(
        summarise_by_location(mod_zetas,mod_post, output="zeta") %>%
          dplyr::rename(!!zeta_label := "zeta"),
        by=c("location")
      ) %>%
      dplyr::inner_join(
        summarise_by_location(mod_times,
                              tibble::tibble(critical_time=c(NA)),
                              output = "critical_time") %>%
          dplyr::rename(!!critical_time_label := "critical_time"),
        by=c("location")
      )

    res_list[[i]] <- res

    i <- i + 1
  }


  result <- res_list %>% purrr::reduce(dplyr::inner_join, by=c("location"))


  return(result)
}
