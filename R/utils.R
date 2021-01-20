

#' Create summary and 95% cri for variable
#' @param v vector
#' @param dp int number of decimal places to format to.
#' @importFrom stats quantile median
#' @noRd
summary_m_95cri <- function(v,dp=2){
  m <- median(v,na.rm = TRUE)
  lc <- quantile(v,0.05,na.rm=TRUE)
  uc <- quantile(v,0.95,na.rm=TRUE)
  glue::glue("{round(m,dp)} ({round(lc,dp)} - {round(uc,dp)})")
}
