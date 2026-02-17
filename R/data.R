#' Sample data for canonical correlation analysis from the SAS manual
#'
#' This dataset includes six variables of fitness club data from the SAS manual.
#'
#' @format A data frame with 20 rows and 6 variables:
#' \describe{
#'   \item{Weight}{Weight measurement}
#'   \item{Waist}{Waist measurement}
#'   \item{Pulse}{Pulse measurement}
#'   \item{Chins}{Number of chin-ups}
#'   \item{Situps}{Number of sit-ups}
#'   \item{Jumps}{Number of jumps}
#' }
#'
#' @source
#' \url{https://documentation.sas.com/doc/en/statcdc/14.2/statug/statug_cancorr_example01.htm}
#'
#' @examples
#' data(sas_ex1)
#'
#' \donttest{
#' ## Canonical Correlation Analysis
#' cancorr(X_vars=c("Weight", "Waist", "Pulse"),
#'         Y_vars=c("Chins", "Situps", "Jumps"),
#'         data=sas_ex1)
#' }
"sas_ex1"

#' Sample data for redundancy analysis from the SAS manual
#'
#' This dataset includes seven variables from the SAS manual.
#'
#' @format A matrix with 10 rows and 7 columns:
#' \describe{
#'   \item{y1, y2, y3}{Y variables}
#'   \item{x1, x2, x3, x4}{X variables}
#' }
#'
#' @source
#' \url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_transreg_details23.htm}
#'
#' @examples
#' data(sas_ex2)
#'
#' \donttest{
#' ## Redundancy Analysis
#' rda(X_vars=c("x1", "x2", "x3", "x4"),
#'     Y_vars=c("y1", "y2", "y3"),
#'     data=sas_ex2)
#' }
"sas_ex2"

#' Correlation matrix of a model of motivation
#'
#' This dataset includes a correlation matrix of 12 variables
#' (n=533) of a model of motivation reported by Chittum, Jones, and
#' Carter (2019).
#'
#' @format A list with the following components:
#' \describe{
#'   \item{data}{A 12x12 correlation matrix.}
#'   \item{n}{A sample size (533).}
#' }
#'
#' @source
#' Chittum, J. R., Jones, B. D., & Carter, D. M. (2019). A person-centered
#' investigation of patterns in college students' perceptions of motivation
#' in a course. *Learning and Individual Differences*, **69**,
#' 94-107. \doi{10.1016/j.lindif.2018.11.007}
#'
#' @references
#'   Gu, F., Yung, Y.-F., Cheung, M. W.-L., Joo, B.-K., & Nimon,
#'   K. (2023). Statistical inference in redundancy analysis: A direct
#'   covariance structure modeling approach. *Multivariate Behavioral
#'   Research*, **58(5)**,
#'   877-893. \doi{10.1080/00273171.2022.2141675}
#'
#' @examples
#' \donttest{
#' data(Chittum19)
#'
#' ## Redundancy Analysis
#' rda(X_vars=c("Empowerment", "Usefulness", "Success", "Interest", "Caring"),
#'     Y_vars=c("Final_Exam", "Learning", "Course_Rating", "Instr_Rating",
#'              "Effort", "Cog_Engage", "Cost"),
#'     Cov=Chittum19$data, numObs=Chittum19$n)
#' }
"Chittum19"

#' Raw data used in Nimon, Joo, and Bontrager (2021)
#'
#' This dataset includes the raw data of 13 variables reported by Nimon, Joo, and Bontrager (2021).
#'
#' @format A data frame with 13 variables.
#'
#' @source
#'   Nimon, K., Joo, B.-K. (Brian), & Bontrager, M. (2021). Work cognitions and work intentions: A canonical correlation study. *Human Resource Development International*, **24(1)**, 65-91. \doi{10.1080/13678868.2020.1775038}
#'
#' @references
#'   Gu, F., & Cheung, M. W.-L. (2023). A Model-based approach to multivariate principal component regression: Selection of principal components and standard error estimates for unstandardized regression coefficients. *British Journal of Mathematical and Statistical Psychology*, **76(3)**, 605-622.
#'   \doi{10.1111/bmsp.12301}
#'
#'   Gu, F., Yung, Y.-F., Cheung, M. W.-L., Joo, B.-K., & Nimon, K.
#'   (2023). Statistical inference in redundancy analysis: A direct
#'   covariance structure modeling
#'   approach. *Multivariate Behavioral Research*, **58(5)**, 877-893. \doi{10.1080/00273171.2022.2141675}
#'
#' @examples
#' \donttest{
#' data(Nimon21)
#'
#' ## Redundancy Analysis
#' rda(X_vars=c("AU", "CC", "CL", "CO", "DF", "FB", "GR", "MW"),
#'     Y_vars=c("IDE", "IEE", "IOCB", "IPR", "ITS"),
#'     data=Nimon21)
#'
#' ## Multivariate Principal Component Regression
#' mpcr(X_vars=c("AU", "CC", "CL", "CO", "DF", "FB", "GR", "MW"),
#'      Y_vars=c("IDE", "IEE", "IOCB", "IPR", "ITS"),
#'      pca="COR", pc_select=1,
#'      data=Nimon21)
#' }
"Nimon21"

#' Correlation matrix of artificial data
#'
#' This dataset includes a correlation matrix of nine artificial
#' variables used in Table 1 of Lambert, Wildt, and Durand (1988).
#'
#' @format A 9x9 correlation matrix.
#'
#' @source
#' Lambert, Z. V., Wildt, A. R., & Durand, R. M. (1988). Redundancy
#' analysis: An alternative to canonical correlation and multivariate
#' multiple regression in exploring interset associations.
#' *Psychological Bulletin*, **104(2)**, 282-289.
#' \doi{10.1037/0033-2909.104.2.282}
#'
#' @references
#' Gu, F., Yung, Y.-F., Cheung, M. W.-L., Joo, B.-K., & Nimon,
#' K. (2023). Statistical inference in redundancy analysis: A direct
#' covariance structure modeling approach. *Multivariate Behavioral
#' Research*,
#' **58(5)**, 877-893. \doi{10.1080/00273171.2022.2141675}
#'
#' @examples
#' data(Lambert88)
#'
#' \donttest{
#' ## Redundancy Analysis
#' rda(X_vars=paste0("x", 1:5), Y_vars=paste0("y", 1:4), Cov=Lambert88, numObs=100)
#' }
"Lambert88"

#' Correlation matrix of a model of disgust
#'
#' This dataset includes a correlation matrix of 13 variables
#' (n=679) between five subscales (y1 to y5) of the Disgust Emotion
#' Scale and eight subscales (x1 to x8) of the Disgust Scale reported by
#' Thorndike (2000, p. 238).
#'
#' @format A list with the following components:
#' \describe{
#'   \item{data}{A 13x13 correlation matrix.}
#'   \item{n}{A sample size (679).}
#' }
#'
#' @source
#' Thorndike, R. M. (2000). Canonical correlation analysis. In
#' H. E. A. Tinsley & S. D. Brown (Eds.), *Handbook of applied multivariate statistics and mathematical modeling*
#' (pp. 237-263). San Diego, CA: Academic Press.
#'
#' @references
#' Gu, F., Yung, Y.-F., & Cheung, M. W.-L. (2019). Four covariance structure models for canonical correlation analysis: A COSAN modeling approach. *Multivariate Behavioral Research*, **54(2)**, 192-223. \doi{10.1080/00273171.2018.1512847}
#'
#' @examples
#' \donttest{
#' data(Thorndike00)
#'
#' ## Canonical Correlation Analysis
#' ## Note. We swap the X_vars and Y_vars because cancorr() expects that
#' ## X_vars cannot have more variables than Y_vars.
#'
#' cancorr(X_vars=c("y1", "y2", "y3", "y4", "y5"),
#'         Y_vars=c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"),
#'         Cov=Thorndike00$data, numObs=Thorndike00$n)
#' }
"Thorndike00"
