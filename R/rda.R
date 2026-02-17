#' Redundancy analysis
#'
#' This function conducts redundancy analysis using the OpenMx package.
#' Missing data are handled with the full information maximum likelihood
#' method when raw data are available. It provides standard errors for
#' the standardized estimates.
#'
#' @param X_vars A vector of characters of the X variables.
#' @param Y_vars A vector of characters of the Y variables.
#' @param data A data frame containing raw data. If NULL, `Cov` and
#'   `numObs` must be provided.
#' @param Cov A covariance or correlation matrix. Required when `data` is NULL.
#' @param numObs A sample size. Required when `data` is NULL.
#' @param extraTries This function calls [OpenMx::mxTryHard()] to obtain
#'   parameter estimates and their standard errors. `extraTries` is the
#'   number of extra runs. If `extraTries=0`, [OpenMx::mxRun()] is
#'   called.
#' @param ... Additional arguments passed to either
#'   [OpenMx::mxTryHard()] or [OpenMx::mxRun()].
#'
#' @return A list with class `RDA`. It stores the model in OpenMx
#'   objects. The fitted object is stored in `mx.fit`.
#'
#' @author Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
#'
#' @references
#'   Gu, F., Yung, Y.-F., Cheung, M. W.-L., Joo, B.-K., & Nimon, K.
#'   (2023). Statistical inference in redundancy analysis: A direct
#'   covariance structure modeling
#'   approach. *Multivariate Behavioral Research*, **58(5)**, 877-893. \doi{10.1080/00273171.2022.2141675}
#'
#' @seealso [Chittum19], [sas_ex2]
#'
#' @export
#'
#' @examples
#' \donttest{
#' ## Redundancy Analysis
#' rda(X_vars=c("x1", "x2", "x3", "x4"),
#'     Y_vars=c("y1", "y2", "y3"),
#'     data=sas_ex2)
#' }
rda <- function(X_vars, Y_vars, data=NULL, Cov=NULL, numObs=NULL, extraTries=50, ...) {
 
    ## Input validation
    if (is.null(data)) {
        if (is.null(Cov)) stop("Either 'data' or 'Cov' must be provided.\n")
        if (is.null(numObs)) stop("'numObs' must be provided when using 'Cov'.\n")
    }
    
    ## Sample Cov as starting values if raw data are given
    if (!is.null(data)) {
        Cov <- cov(data, use="pairwise.complete.obs")
    }
  
    ## No. of X variables
    p <- length(X_vars)

    ## No. of Y variables
    q <- length(Y_vars)

    ## Get the R matrices from the data
    R <- cov2cor(Cov)
    ## Keep matrix shape for single-variable Y (q==1); vector dropping
    ## causes non-conformable matrix multiplications below.
    Rxx <- R[X_vars, X_vars, drop=FALSE]
    Ryy <- R[Y_vars, Y_vars, drop=FALSE]
    Rxy <- R[Y_vars, X_vars, drop=FALSE]

    tri_prod  <- solve(Rxx) %*% t(Rxy) %*% Rxy
    ei <- eigen(tri_prod)
    ## Take the real part of the eigenvectors (if there are complex)
    vec <- Re(ei$vectors)
    
    ## Check if d=p-q>=2
    d <- p-q
    
    ## D_s: A pxp matrix of T/F elements.
    ## The last d_s=d*(d-1)/2 elements are arbitrarily fixed at 0
    if (d>=2) {
      d_s <- d*(d-1)/2
      D_s <- matrix(rep(c(TRUE, FALSE), times=c(p^2-d_s, d_s)), nrow=p, ncol=p)
      
      ## Don't fix the last d_s elements to 0
      ## Use the original starting values
      ## ei$vectors[!D_s] <- 0
    } else {
      ## All elements are free
      D_s <- matrix(TRUE, nrow=p, ncol=p)
    }

    ## scale is a diagonal matrix
    ## scale <- diag(sqrt(diag(t(vec) %*% Rxx %*% vec)))
    ## start_W <- vec %*% solve(scale)

    ## Starting values for W    
    scale <- sqrt(diag(t(vec) %*% Rxx %*% vec))
    start_W <- vec %*% diag(1/(scale))   

    ## Starting values for Ly
    start_Ly <- Rxy %*% start_W  
  
    ## Standard deviations of X and Y
    SD <- mxMatrix("Diag", nrow=(p+q), ncol=(p+q), free=TRUE, 
                   values=sqrt(diag(Cov[c(X_vars, Y_vars), c(X_vars, Y_vars), drop=FALSE])), 
                   labels=c(paste0("sdx", 1:p), paste0("sdy", 1:q)),
                   name="SD")
  
    ## qxp matrix of zero
    qZerop <- mxMatrix("Zero", nrow=q, ncol=p, name="qZerop")
  
    ## Identity matrices in dimensions of p and q
    Iq <- mxMatrix("Iden", nrow=q, ncol=q, name="Iq") 
    Ip <- mxMatrix("Iden", nrow=p, ncol=p, name="Ip") 

    ## W matrix with random starting values
    W <- mxMatrix("Full", nrow=p, ncol=p, free=D_s,
                  values=start_W,
                  labels=outer(1:p, 1:p, function(x, y) paste0("W", x, "_", y)),
                  name="W")

    ## Block diagonal matrix of W and I
    W_I <- mxAlgebra(rbind(cbind(solve(t(W)), t(qZerop)),
                           cbind(qZerop, Iq)), name="W_I")

    ## Ryy
    ## Keep 1x1 matrix form when q==1 so cov2cor()/vechs() remain valid.
    Ryy <- mxMatrix("Stand", nrow=q, ncol=q, free=TRUE,
                    values=vechs(cov2cor(Cov[Y_vars, Y_vars, drop=FALSE])),
                    labels=vechs(outer(1:q, 1:q, function(x, y) paste0("Ry", x, "_", y))),
                    name="Ryy")

    ## Ly
    if (d > 0) {
        ## Ly when p > q, the last (p-q) columns are zeros
        Ly <- mxMatrix("Full", nrow=q, ncol=p,
                       free=cbind(matrix(TRUE, nrow=q, ncol=q), 
                                  matrix(FALSE, nrow=q, ncol=(p-q))),
                       values=cbind(matrix(start_Ly[1:q, 1:q], nrow=q, ncol=q), 
                                    matrix(0, nrow=q, ncol=(p-q))),
                       labels=outer(1:q, 1:p, function(x, y) paste0("Ly", x, "_", y)),
                       name="Ly")
    } else {
        ## Ly when p <= q
        Ly <- mxMatrix("Full", nrow=q, ncol=p,
                       free=TRUE,
                       values=start_Ly,
                       labels=outer(1:q, 1:p, function(x, y) paste0("Ly", x, "_", y)),
                       name="Ly")
    }

    ## Matrix to store q, the no. of Y variables
    q_matrix <- mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=q, name="q")

    ## Lambdas: Equation (7)
    Lambdas <- mxAlgebra(diag2vec(t(Ly) %*% Ly)/q, name="Lambdas")

    ## Matrix for cumulative sum
    K <- matrix(1, nrow=p, ncol=p)
    K[upper.tri(K)] <- 0
    K <- mxMatrix("Full", nrow=p, ncol=p, free=FALSE, values=K, name="K")
    
    ## Cumulative sum of Lambdas
    Lambdas_cum <- mxAlgebra(K %*% Lambdas, name="Lambdas_cum")
    
    ## Lx: Equation (12)
    Lx <- mxAlgebra(solve(t(W)), name="Lx")

    ## Block matrix of I and R
    I_R <- mxAlgebra(rbind(cbind(Ip, t(Ly)),
                           cbind(Ly, Ryy)), name="I_R")

    ## px1 column vector of ones
    pOne1 <- mxMatrix("Unit", nrow=p, ncol=1, name="pOne1")
  
    ## p*(p-1)/2x1 column vector of zeros 
    psZero1 <- mxMatrix("Zero", nrow=p*(p-1)/2, ncol=1, name="psZero1")

    ## Expected covariance matrix: Equation (10)
    expCov <- mxAlgebra(SD %&% (W_I %&% I_R), name="expCov")

    ## Equation (11)
    constraint1 <- mxConstraint( diag2vec(solve(t(W)) %*% solve(W)) == pOne1, 
                                name = 'constraint1' )

    ## Second constraint that off diagonals are zeros 
    constraint2 <- mxConstraint( vechs(t(Ly) %*% Ly) == psZero1, name = 'constraint2' )

    if (is.null(data)) {
        ## Cov or Cor as inputs
        mxdata <- mxData(observed=Cov[c(X_vars, Y_vars), c(X_vars, Y_vars), drop=FALSE], 
                         type="cov", numObs=numObs)
    
        ## Expected covariance, which is required with raw data, in the fit function
        expFun <- mxExpectationNormal(covariance="expCov", dimnames=c(X_vars, Y_vars))
    
        ## Combine everything to a model
        mx.model <- mxModel("RA", mxdata,
                            qZerop, SD, Iq, Ip, W, W_I, Ryy, Ly, I_R, pOne1, psZero1, Lambdas, 
                            K, Lambdas_cum, Lx, q_matrix, expCov, expFun,
                            mxFitFunctionML(), constraint1, constraint2)
  
    } else {
      ## Raw data as inputs
        mxdata <- mxData(observed=data[, c(X_vars, Y_vars)], type="raw")
    
        expMean <- mxMatrix("Full", nrow=1, ncol=p+q, free=TRUE, 
                            values=apply(data[, c(X_vars, Y_vars)], 2, mean, na.rm=TRUE),
                            labels = c(paste0("mux_", seq_len(p)), paste0("muy_", seq_len(q))),
                            name="expMean")
    
        ## Expected covariance in the fit function
        expFun <- mxExpectationNormal(covariance="expCov", means="expMean", 
                                      dimnames=c(X_vars, Y_vars))    
    
        ## Combine everything to a model
        mx.model <- mxModel("RA", mxdata,
                            qZerop, SD, Iq, Ip, W, W_I, Ryy, Ly, I_R, pOne1, psZero1, Lambdas,
                            K, Lambdas_cum, Lx, q_matrix, expCov, expFun, expMean, 
                            mxFitFunctionML(), constraint1, constraint2)
    }

    ## Use the starting values as the final estimates.
    ## Do not activate the optimizer
    ## if (optimizer==FALSE) {
    ##     plan <- omxDefaultComputePlan()
    ##     plan$steps <- list(plan$steps$ND, plan$steps$SE, plan$steps$RD, plan$steps$RE)
    ##     mx.model <- mxModel(mx.model, plan)
    ## }

    if (extraTries==0) {
        mx.fit <- suppressMessages(mxRun(mx.model, ...))
    } else {
        mx.fit <- mxTryHard(mx.model, extraTries = extraTries, ...) 
        ## Run it one more time to minimize error
        mx.fit <- suppressMessages(mxRun(mx.fit))
    }
    
    #### Constraints for checking
    ## Diagonals should be 1: Equation (14)
    Constraint1 <- mxEval(t(diag(Lx %*% t(Lx))), mx.fit)
    colnames(Constraint1) <- X_vars

    ## Off-diagonals should be 0
    Constraint2 <- mxEval(vechs(t(Ly) %*% Ly), mx.fit)
    Constraint2 <- c(min(Constraint2), max(Constraint2))
    names(Constraint2) <- c("Min", "Max")

    ## W matrix
    W_est <- mxEval(W, mx.fit)
    dimnames(W_est) <- list(X_vars, X_vars)

    ## Ly matrix
    Ly_est <- mxEval(Ly, mx.fit)
    dimnames(Ly_est) <- list(Y_vars, X_vars)

    ## Lx matrix
    Lx_est <- mxEval(Lx, mx.fit)
    dimnames(Lx_est) <- list(X_vars, X_vars)

    ## Lambdas
    Lambdas_est <- mxEval(Lambdas, mx.fit)
    Lambdas_SE <- suppressMessages(mxSE(Lambdas, mx.fit))

    ## Lambdas cumulative sum
    Lambdas_cum_est <- mxEval(Lambdas_cum, mx.fit)
    Lambdas_cum_SE <- suppressMessages(mxSE(Lambdas_cum, mx.fit))  
    
    ## SE of W matrix
    W_SE <- suppressMessages(mxSE(W, mx.fit))
    dimnames(W_SE) <- list(X_vars, X_vars)
  
    ## SE of Ly matrix
    Ly_SE <- suppressMessages(mxSE(Ly, mx.fit))
    dimnames(Ly_SE) <- list(Y_vars, X_vars)    
  
    ## SE of Lx matrix
    Lx_SE <- suppressMessages(mxSE(Lx, mx.fit))
    dimnames(Lx_SE) <- list(X_vars, X_vars)   
    
    out <- list(mx.fit=mx.fit, mx.model=mx.model, Constraint1=Constraint1,
                Constraint2=Constraint2, W_est=W_est, Ly_est=Ly_est,
                Lx_est=Lx_est, Lambdas_est=Lambdas_est, Lambdas_cum_est=Lambdas_cum_est, 
                W_SE=W_SE, Ly_SE=Ly_SE, Lx_SE=Lx_SE, Lambdas_SE=Lambdas_SE,
                Lambdas_cum_SE=Lambdas_cum_SE)  
    class(out) <- "RDA"
  
    # options(warn = oldw)
    out
}

#' Print Method for RDA Objects
#'
#' Print method for `RDA` objects.
#'
#' @param x An object returned from the class of `RDA`.
#' @param digits Number of digits in printing the matrices. The default is 4.
#' @param ... Unused.
#'
#' @return No return value, called for side effects
#'
#' @author Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
#'
#' @method print RDA
#' @export
print.RDA <- function(x, digits=4, ...) {
    if (!is.element("RDA", class(x)))
        stop("\"x\" must be an object of class \"RDA\".")

    cat("Please check the constraints before interpreting the results.\n")
    cat("Constraint 1. The following values should be close to 1:", round(x$Constraint1, 6), ".\n")
    cat("Constraint 2. The following values (min and max) should be close to 0:", round(x$Constraint2, 6), ".\n")
  
    cat("\nW matrix:\n")
    .mprint(x$W_est, digits=digits)
    cat("\nW matrix (SE):\n")
    .mprint(x$W_SE, digits=digits)
  
    cat("\nLy matrix:\n")
    .mprint(x$Ly_est, digits=digits)
    cat("\nLy matrix (SE):\n")
    .mprint(x$Ly_SE, digits=digits)
  
    cat("\nLx matrix:\n")
    .mprint(x$Lx_est, digits=digits)
    cat("\nLx matrix (SE):\n")
    .mprint(x$Lx_SE, digits=digits)
  
    cat("\nIndividual redundancy:\n")
    .mprint(x$Lambdas_est, digits=digits)
    cat("\nIndividual redundancy (SE):\n")
    .mprint(x$Lambdas_SE, digits=digits)

    cat("\nCumulative redundancy:\n")
    .mprint(x$Lambdas_cum_est, digits=digits)
    cat("\nCumulative redundancy (SE):\n")
    .mprint(x$Lambdas_cum_SE, digits=digits)
}
