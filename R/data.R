#' Clinical trial data from 369 Alzheimer's disease patients.
#' 
#' A dataset containing the treatment assignment, baseline measurements,
#' and outcomes of 369 patients enrolled in four clinical trials
#' for treatments of Alzheimer's disease. Patients assigned to the test
#' treatments are excluded, leaving patients assigned to either a placebo
#' or the palliative standard of care.
#' 
#' @format A data frame with 369 rows and 6 variables:
#' \describe{
#'   \item{Treatment}{treatment assignment}
#'   \item{Severity}{baseline severity (ADAS-Cog 11)}
#'   \item{Decline}{decline in score on Mini Mental State Exam (MMSE) divided
#'                  by number of years since onset of symptoms}
#'   \item{Sex}{male or female}
#'   \item{Carrier}{carrier status of ApoE4 allele}
#'   \item{Improvement}{baseline severity minus severity at 12 weeks
#'                      (positive is good)}
#' }
#' 
#' @source AbbVie, Inc.
"alzheimers"