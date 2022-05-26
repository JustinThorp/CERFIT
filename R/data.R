#' Observational Educational Dataset
#'
#' A simulated dataset containing the grades and other attributes of 1000
#' simulated students
#'
#'
#' @format A data frame with 1000 rows and 7 variables:
#' \describe{
#'   \item{SAT_MATH}{SAT Math Score}
#'   \item{HSGPA}{High School GPA}
#'   \item{AGE}{Age of Student}
#'   \item{GENDER}{Gender of Student}
#'   \item{URM}{Under Represented Minority}
#'   \item{A}{Treatment Variable}
#'   \item{Y}{Students Final Grade}
#' }
#' @source Wilke, Morten C., et al. “Estimating the Optimal Treatment
#' Regime for Student Success Programs.” Behaviormetrika, vol. 48, no. 2, 2021,
#' pp. 309–343., https://doi.org/10.1007/s41237-021-00140-0.
"educational"


#' Randomized Controlled Trial Warts Dataset
#'
#' A dataset comparing immunotherapy to cryotherapy treatments and their effeteness of
#' removing warts
#'
#'
#' @format A data frame with 180 rows and 8 variables:
#' \describe{
#'   \item{sex}{Patients Sex}
#'   \item{age}{Patients Age}
#'   \item{Time}{Time Elapsed Before Treatment}
#'   \item{Number_of_Warts}{Number of Warts}
#'   \item{Type}{Type of Wart}
#'   \item{Area}{Wart Surface Area}
#'   \item{Result_of_Treatment}{Treatment Outcome}
#'   \item{treatment}{0 for immunotherapy and 1 for cryotherapy}
#' }
#' @source Khozeimeh, Fahime, et al. “An Expert System for Selecting Wart
#' Treatment Method.” Computers in Biology and Medicine, vol. 81, 2017,
#' pp. 167–175., https://doi.org/10.1016/j.compbiomed.2017.01.001.
"warts"
