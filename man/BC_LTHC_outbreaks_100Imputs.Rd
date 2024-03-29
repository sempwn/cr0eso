% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BC_LTHC_outbreaks_100Imputs}
\alias{BC_LTHC_outbreaks_100Imputs}
\title{Long Term Health Care facility COVID-19 outbreaks in British Columbia, with 100 imputations of missing data}
\format{
A list with 100 elements, corresponding to 100 imputations of missing symptom onset data, each with 6 sub-elements
\describe{
  \item{Location}{Long Term Health Care facility identifier}
  \item{num_cases}{Total number of COVID-19 cases identified in each location's outbreak}
  \item{time_series}{Daily symptom onset incidence in each location's outbreak e.g. time_series[[2]][3] gives the number of cases with symptom onset on the third day of Location[2]'s outbreak.}
  \item{capacity}{Resident capacity of each location}
  \item{reported_date}{Date of first case reported in each location's outbreak}
  \item{case_matrix}{time_series objects combined into a matrix, where each column is a different facility outbreak and rows 1,...,n give the number of cases with symptom onset on day 1,...,n in that outbreak. Used in Stan model fitting.}
}
}
\source{
British Columbia Centre for Disease Control
}
\usage{
BC_LTHC_outbreaks_100Imputs
}
\description{
A dataset containing 100 replicates of LTHC covid-19 outbreaks in BC, Canada, each
with a random imputation of missing symptom onset times.
}
\keyword{datasets}
