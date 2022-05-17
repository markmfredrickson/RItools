#' The ASSIST Trial baseline data from Yudkin and Moher 2001 consist of 21
#' general practices containing 2142 patients used for the design of a
#' randomized trial which assigned to three treatments aiming to compare methods
#' of preventing coronary heart disease. We have expanded the aggregated data
#' from the practice level to the individual level and added a simulated
#' randomized treatment variable.
#'
#' @source The data come from Table II on page 345 of Yudkin and Moher 2001,
#' Statistics in Medicine.
#'
#' @title ASSIST Trial Data from Yudkin and Moher 2001
#'
#' @format A data frame with  2142 rows and 9 columns
#'
#' \itemize{
#'   \item practice: Identifier for the practice
#'   \item id: Identifier for patients
#'   \item n_practice: Number of patients with CHD in the practice
#'   \item assessed: Whether a patient was coded as "adequately assessed" (the
#'        outcome of the study, measured here at baseline).
#'   \item aspirin: Whether a patient was treated with aspirin at baseline.
#'   \item hypo: Whether a patient was treated with hypotensives at
#'        baseline.
#'   \item lipid: Whether patient was treated with lipid-lowering drugs
#'        at baseline.
#'   \item assess_strata: Strata of the practice defined by the proportion of
#'        people adequately assessed at baseline (three strata, following Yudkin
#'        and Moher 2001).
#'   \item trt: A simulated binary treatment, assigned at random at the practice
#'        level within levels of assess_strata.
#'}
#' @references Yudkin, P. L. and Moher, M. 2001. "Putting theory into practice:
#' a cluster randomized trial with a small number of clusters" \emph{Statistics
#' in Medicine}, 20:341-349.
#' @keywords datasets
"ym_long"



