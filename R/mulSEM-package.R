#' @keywords internal
"_PACKAGE"

#' mulSEM: Some Multivariate Analyses using Structural Equation Modeling
#'
#' A set of functions for some multivariate analyses utilizing a
#' structural equation modeling (SEM) approach through the 'OpenMx' package.
#' These analyses include canonical correlation analysis (CANCORR),
#' redundancy analysis (RDA), and multivariate principal component
#' regression (MPCR). It implements procedures discussed in
#' Gu and Cheung (2023) <doi:10.1111/bmsp.12301>, Gu, Yung, and Cheung
#' (2019) <doi:10.1080/00273171.2018.1512847>, and Gu et al. (2023)
#' <doi:10.1080/00273171.2022.2141675>.
#'
#' @author Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>, Fei Gu <fgu@vt.edu>,
#'   Yiu-Fai Yung <Yiu-Fai.Yung@sas.com>
#'
#' @references
#' Gu, F., & Cheung, M. W.-L. (2023). A model-based approach to
#' multivariate principal component regression: Selection of principal
#' components and standard error estimates for unstandardized regression
#' coefficients. *British Journal of Mathematical and Statistical
#' Psychology*, **76(3)**, 605-622. \doi{10.1111/bmsp.12301}
#'
#' Gu, F., Yung, Y.-F., & Cheung, M. W.-L. (2019). Four covariance
#' structure models for canonical correlation analysis: A COSAN modeling
#' approach. *Multivariate Behavioral Research*, **54(2)**,
#' 192-223. \doi{10.1080/00273171.2018.1512847}
#'
#' Gu, F., Yung, Y.-F., Cheung, M. W.-L., Joo, B.-K., & Nimon,
#' K. (2023). Statistical inference in redundancy analysis: A direct
#' covariance structure modeling approach. *Multivariate Behavioral
#' Research*, **58(5)**, 877-893. \doi{10.1080/00273171.2022.2141675}
#'
#' @examples
#' \donttest{
#' ## Canonical Correlation Analysis
#' cancorr(X_vars=c("Weight", "Waist", "Pulse"),
#'         Y_vars=c("Chins", "Situps", "Jumps"),
#'         data=sas_ex1)
#'
#' ## Redundancy Analysis
#' rda(X_vars=c("x1", "x2", "x3", "x4"),
#'     Y_vars=c("y1", "y2", "y3"),
#'     data=sas_ex2)
#'
#' ## Multivariate Principal Component Regression
#' mpcr(X_vars=c("AU", "CC", "CL", "CO", "DF", "FB", "GR", "MW"),
#'      Y_vars=c("IDE", "IEE", "IOCB", "IPR", "ITS"),
#'      pca="COR", pc_select=1,
#'      data=Nimon21)
#' }
#'
#' @name mulSEM-package
#' @aliases mulSEM
NULL
