#' Long Term Health Care facility COVID-19 outbreaks in British Columbia, with 100 imputations of missing data
#'
#' A dataset containing 100 replicates of LTHC covid-19 outbreaks in BC, Canada, each
#' with a random imputation of missing symptom onset times.
#'
#'
#' @format A list with 100 elements, corresponding to 100 imputations of missing symptom onset data, each with 6 sub-elements
#' \describe{
#'   \item{Location}{Long Term Health Care facility identifier}
#'   \item{num_cases}{Total number of COVID-19 cases identified in each location's outbreak}
#'   \item{time_series}{Daily symptom onset incidence in each location's outbreak e.g. time_series[[2]][3] gives the number of cases with symptom onset on the third day of Location[2]'s outbreak.}
#'   \item{capacity}{Resident capacity of each location}
#'   \item{reported_date}{Date of first case reported in each location's outbreak}
#'   \item{case_matrix}{time_series objects combined into a matrix, where each column is a different facility outbreak and rows 1,...,n give the number of cases with symptom onset on day 1,...,n in that outbreak. Used in Stan model fitting.}
#'   ...
#' }
#' @source British Columbia Centre for Disease Control
"BC_LTHC_outbreaks_100Imputs"



#' Additional covariate data on Long Term Health Care facility COVID-19 outbreaks in British Columbia
#'
#' A dataset containing additional covariate data for each of the 18 Long Term Health Care
#' facility COVID-19 outbreaks in British Columbia contained in dataset
#' 'BC_LTHC_outbreaks_100Imputs'. The majority of this data was sourced from the Office of
#' the Seniors Advocate, British Columbia, Long-Term Care Facilities Quick Facts Directory.
#' https://www.seniorsadvocatebc.ca/quickfacts/location. All OSABC data is from the 2018/19
#' year. Identity of the initial case in each facility outbreak was obtained from the BCCDC.
#'
#'
#' @format A dataframe with 18 rows, corresponding to the 18 LTHC outbreaks in this dataset, and 12 columns of data concerning each outbreak
#' \describe{
#'   \item{Location}{Long Term Health Care facility identifier}
#'   \item{LTHC facility type}{Funding type of each facility: Private, Private non-profit or Public, as much as was possible to conclude}
#'   \item{Year facility opened}{Year in which each facility first opened to residents}
#'   \item{Resident room type}{Type of accomodations available for residents: single rooms only, multi-resident rooms only, or a mix (either primarily single rooms or primarily multi-person rooms)}
#'   \item{Accreditation status}{Has each facility been voluntarily accredited by Accreditation Canada (see https://www2.gov.bc.ca/gov/content/health/accessing-health-care/home-community-care/accountability/quality-and-safety)}
#'   \item{Direct care hours /resident/day}{Total direct care hours (hours per resident per day), Nursing/Care Aide care + Allied health hours, in each facility}
#'   \item{Average resident age (years)}{Average age in years of residents in each facility}
#'   \item{Average resident stay (days)}{Average length of stay in days of residents in each facility}
#'   \item{Residents dependent for daily activities (%)}{Percent of residents in each facility who are totally dependent in their in activities of daily living}
#'   \item{Number of lodged complaints 2018/19}{Number of licensing complaints lodged in each facility during the 2018/19 year.}
#'   \item{Number of disease outbreaks 2018/19}{Total number of disease outbreak or occurrence incidents recorded in each facility during the 2018/19 year}
#'   \item{Identity of initial COVID-19 case}{Identity (staff/worker or resident/patient) of the initial COVID-19 case in each facility, as determined by earliest recorded date of symptom onset.}
#'   ...
#' }
#' @source Office of the Seniors Advocate, British Columbia \url{https://www.seniorsadvocatebc.ca/quickfacts/location} and British Columbia Centre for Disease Control
"BC_OSABC_facilitydata"

#' Object created from `seir_model_fit` for BC data with fixed intervention. See vignettes for more information.
"bc_fit"

#' Object created from `seir_model_fit` for BC data with hierarchical intervention. See vignettes for more intervention.
"bc_fit_zeta"
