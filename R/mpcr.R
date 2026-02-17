#' Multivariate Principal Component Regression (MPCR)
#'
#' Conduct multivariate principal component regression
#'
#' This function conducts multivariate principal component regression
#' using the OpenMx package. Missing data are handled with the full
#' information maximum likelihood method when raw data are available.
#' It provides standard errors for the estimates.
#'
#' @param X_vars A vector of characters of the X variables.
#' @param Y_vars A vector of characters of the Y variables.
#' @param data A data frame containing raw data. If NULL, `Cov` and
#'   `numObs` must be provided.
#' @param Cov A covariance or correlation matrix. Required when `data` is NULL.
#' @param Means An optional mean vector. Can be provided when `data` is NULL.
#' @param numObs A sample size. Required when `data` is NULL.
#' @param pca Whether principal component analysis is based on
#'   unstandardized (`COV`) or standardized (`COR`) variables.
#' @param pc_select PCs selected in the regression analysis. For example,
#'   `pc_select=c(1,2)` to use the first two PCs in the multiple
#'   regression analysis.
#' @param extraTries This function calls [OpenMx::mxTryHard()] to obtain
#'   parameter estimates and their standard errors. `extraTries` is the
#'   number of extra runs. If `extraTries=0`, [OpenMx::mxRun()] is
#'   called.
#' @param ... Additional arguments passed to either
#'   [OpenMx::mxTryHard()] or [OpenMx::mxRun()].
#'
#' @return A list with class `MPCR`. It stores the model in OpenMx
#'   objects. The fitted object is stored in `mx.fit`.
#'
#' @author Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
#'
#' @references
#'   Gu, F., & Cheung, M. W.-L. (2023). A model-based approach to
#'   multivariate principal component regression: Selection of principal
#'   components and standard error estimates for unstandardized regression
#'   coefficients. *British Journal of Mathematical and Statistical
#'   Psychology*, **76(3)**, 605-622.
#'   \doi{10.1111/bmsp.12301}
#'
#' @seealso [Nimon21]
#'
#' @export
#'
#' @examples
#' \donttest{
#' ## Multivariate Principal Component Regression
#' mpcr(X_vars=c("AU", "CC", "CL", "CO", "DF", "FB", "GR", "MW"),
#'      Y_vars=c("IDE", "IEE", "IOCB", "IPR", "ITS"),
#'      pca="COR", pc_select=1,
#'      data=Nimon21)
#' }
mpcr <- function(X_vars, Y_vars, data=NULL, Cov=NULL, Means=NULL, numObs=NULL, pca=c("COV", "COR"),
                 pc_select=NULL, extraTries=50, ...) {

    pca <- match.arg(pca)
    
    ## Input validation
    if (is.null(data)) {
        if (is.null(Cov)) stop("Either 'data' or 'Cov' must be provided.\n")
        if (is.null(numObs)) stop("'numObs' must be provided when using 'Cov'.\n")
    }
    
    ## Some simple check on pc_select
    if (is.null(pc_select)) stop("Please indicate which PCs you want to use in the multiple regression analysis, e.g., 'pc_select=c(1,2)' to use the first two PCs.\n")    
    pc_select <- sort(pc_select)
    if (!all(pc_select %in% seq_along(X_vars))) stop("The PC selected must be between 1 and ", length(X_vars), ".\n") 
    if (length(pc_select) != length(unique(pc_select))) stop("Some of the PC are repeated.\n")
    
    ## Whether the means are given
    if (!is.null(data) || !is.null(Means)) {
        mean.structure <- TRUE
    } else {
        mean.structure <- FALSE
    }

    ## Sample Cov as starting values if raw data are given
    if (!is.null(data)) {
        data <- data[, c(X_vars, Y_vars)]
        Cov <- cov(data, use="pairwise.complete.obs")

        ## Raw data as inputs
        mxdata <- mxData(observed=data, type="raw")
    } else {
        ## Cov or Cor as inputs
        mxdata <- mxData(observed=Cov[c(X_vars, Y_vars), c(X_vars, Y_vars)], 
                         type="cov", numObs=numObs)
    }
  
    ## No. of X variables
    p <- length(X_vars)

    ## No. of Y variables
    q <- length(Y_vars)

    ## Starting values of eigen values decomposition on either COV or COR
    if (pca=="COV") {
        S <- Cov        
    } else {
        ## Correlation structure
        S <- cov2cor(Cov)
    }
    ## Keep matrix shape when p==1; eigen() on a dropped scalar is invalid.
    eig <- eigen(S[X_vars, X_vars, drop=FALSE])
    ## V_start <- eig$vectors 
    ## Lambda_start <- diag(eig$values)
    ## Keep matrix shapes for q==1/p==1 to preserve conformable products.
    H_start <- S[Y_vars, X_vars, drop=FALSE] %*% solve(t(eig$vectors))

    ## Prepare some matrices
    Lambda <- mxMatrix(type="Diag", ncol=p, nrow=p, free=TRUE, values=eig$values,
                       labels=paste0("lambda", seq_len(p)), name="Lambda")

    V <- mxMatrix(type="Full", nrow=p, ncol=p, free=TRUE, values=eig$vectors,
                  labels=outer(seq_len(p), seq_len(p), function(x, y) paste0("v", x, "_", y)),
                  name="V")

    ## A pxp identity matrix
    Ip <- mxMatrix(type="Iden", nrow=p, ncol=p, name="Ip")

    ## A 1xp matrix of ones
    Oneq <- mxMatrix(type="Unit", nrow=q, ncol=1, name="Oneq")
     
    ## Constraint on V
    constraint1 <- mxConstraint(vech(V%*%t(V)) == vech(Ip), name="constraint1")

    ## H is q x p, so labels must also be generated in q x p order.
    ## Using p x q indices here would mismatch dimensions.
    H <- mxMatrix(type="Full", nrow=q, ncol=p, free=TRUE, values=H_start,
                  labels=outer(seq_len(q), seq_len(p), function(x, y) paste0("h", x, "_", y)),
                  name="H")

    if (pca=="COV") {
        ## Keep 1x1 covariance matrix form when q==1.
        Psi <- mxMatrix(type="Symm", nrow=q, ncol=q, free=TRUE, values=vech(Cov[Y_vars, Y_vars, drop=FALSE]),
                        labels=vech(outer(seq_len(q), seq_len(q), function(x, y) paste0("psi", x, "_", y))),
                        name="Psi")
    } else {
        ## Keep 1x1 correlation matrix form when q==1.
        Psi <- mxMatrix(type="Stand", nrow=q, ncol=q, free=TRUE, values=vechs(cov2cor(Cov[Y_vars, Y_vars, drop=FALSE])),
                        labels=vechs(outer(seq_len(q), seq_len(q), function(x, y) paste0("psi", x, "_", y))),
                        name="Psi")
    }        

    ## Mean structure
    if (mean.structure) {

        ## Mean structure if with raw data or Means are provided
        if (!is.null(data)) {
            ## drop=FALSE keeps apply() in matrix mode for single X/Y variables.
            X_means <- matrix(apply(data[, X_vars, drop=FALSE], 2, mean, na.rm=TRUE), ncol=1)
            Tau_start <- apply(data[, Y_vars, drop=FALSE], 2, mean, na.rm=TRUE)
        } else if (!is.null(Means)) {
            X_means <- matrix(Means[X_vars], ncol=1)
            Tau_start <- Means[Y_vars]
        }

        ## Further adjustment if it is correlation structure
        if (pca=="COV") {
            Alpha_start <- solve(eig$vectors) %*% X_means
        } else {            
            Alpha_start <- solve(eig$vectors) %*% (X_means/sqrt(diag(Cov[X_vars, X_vars, drop=FALSE])))
            Tau_start <- Tau_start/sqrt(diag(Cov[Y_vars, Y_vars, drop=FALSE]))
        }
        
        Alpha <- mxMatrix(type="Full", nrow=p, ncol=1, free=TRUE, 
                          values=Alpha_start, labels=paste0("alpha", seq_len(p)), name="Alpha")

        Tau <- mxMatrix(type="Full", nrow=q, ncol=1, free=TRUE, 
                        values=Tau_start, labels=paste0("tau", seq_len(q)), name="Tau")           
    }

    if (pca=="COV") {            
        ## Covariance structure: Equation 3
        expCov <- mxAlgebra( rbind(cbind(V %*% Lambda %*% t(V), V %*% t(H)),
                                   cbind(H %*% t(V), Psi)), name="expCov" )

        ## With mean structure
        if (mean.structure) {

            ## Equation 3
            expMean <- mxAlgebra(cbind(t(V %*% Alpha), t(Tau)), name="expMean")
            
            ## Expected covariance in the fit function
            expFun <- mxExpectationNormal(covariance="expCov", means="expMean", 
                                          dimnames=c(X_vars, Y_vars))    
    
            ## Combine everything to a model
            mx.model <- mxModel("MPCR", mxdata, Lambda, V, H, Psi, Alpha, Tau, Ip, Oneq,
                                constraint1, expCov, expMean, expFun, mxFitFunctionML())
        } else {
            ## No mean structure
            ## Expected covariance in the fit function
            expFun <- mxExpectationNormal(covariance="expCov", dimnames=c(X_vars, Y_vars))    
    
            ## Combine everything to a model
            mx.model <- mxModel("MPCR", mxdata, Lambda, V, H, Psi, Ip, Oneq,
                                constraint1, expCov, expFun, mxFitFunctionML())
        }        
    } else {
        ## Correlation structure
        
        ## Standard deviations of X and Y
        ## Keep 1x1 diagonal matrix form for p==1/q==1.
        Dx <- mxMatrix("Diag", nrow=p, ncol=p, free=TRUE, values=sqrt(diag(Cov[X_vars, X_vars, drop=FALSE])), 
                        labels=paste0("sdx", 1:p),  name="Dx")

        Dy <- mxMatrix("Diag", nrow=q, ncol=q, free=TRUE, values=sqrt(diag(Cov[Y_vars, Y_vars, drop=FALSE])), 
                         labels=paste0("sdy", 1:q),  name="Dy")

        pZeroq <- mxMatrix("Full", nrow=p, ncol=q, free=FALSE, name="pZeroq")

        SD <- mxAlgebra( rbind(cbind(Dx, pZeroq),
                               cbind(t(pZeroq), Dy)), name="SD")
        
        ## SD <- mxMatrix("Diag", nrow=(p+q), ncol=(p+q), free=TRUE, 
        ##                values=sqrt(diag(Cov[c(X_vars, Y_vars), c(X_vars, Y_vars)])), 
        ##                labels=c(paste0("sdx", 1:p), paste0("sdy", 1:q)),  name="SD")

        ## Equation 4
        expCov <- mxAlgebra( SD %&% rbind(cbind(V%*%Lambda%*%t(V), V%*%t(H)),
                                          cbind(H%*%t(V), Psi)), name="expCov" )

        ## Constraint on V %&% Lambda for correlation structure
        constraint2 <- mxConstraint( diag2vec(V %&% Lambda) == diag2vec(Ip), name="constraint2")
        
        ## With mean structure
        if (mean.structure) {

            ## Equation 4
            expMean <- mxAlgebra(t(SD %*% rbind(V %*% Alpha, Tau)), name="expMean")
            
            ## Expected covariance in the fit function
            expFun <- mxExpectationNormal(covariance="expCov", means="expMean",
                                          dimnames=c(X_vars, Y_vars))    
    
            ## Combine everything to a model
            mx.model <- mxModel("MPCR", mxdata, Lambda, V, H, Psi, Alpha, Tau, Ip, Oneq, SD, Dx, Dy,
                                pZeroq, constraint1, constraint2, expCov, expMean, expFun,
                                mxFitFunctionML())
        } else {
            ## No mean structure in this branch: do not include expMean.
            ## Expected covariance in the fit function
            expFun <- mxExpectationNormal(covariance="expCov", dimnames=c(X_vars, Y_vars))

            ## Combine everything to a model
            mx.model <- mxModel("MPCR", mxdata, Lambda, V, H, Psi, Ip, Oneq, SD, Dx, Dy,
                                pZeroq, constraint1, constraint2, expCov, expFun,
                                mxFitFunctionML())
        }
    }

    ## Selection matrix: pxm
    Fpm  <- mxMatrix("Full", nrow=length(X_vars), ncol=length(pc_select), free=FALSE,
                     values=diag(length(X_vars))[, pc_select], name="Fpm")
    
    ## ## Prediction equations
    if (pca=="COV") {
        ## Equation 5
        J <- mxAlgebra( chol(solve(vec2diag(diag2vec(Psi)))) %*% H %*%
                        chol(solve(vec2diag(diag2vec(Lambda)))), name="J")
        ## Eq. 5 uses diag(Psi), not diag(Psi)^{-1}, in the pi weighting term.
        pi <- mxAlgebra( (t(J*J) %*% vec2diag(diag2vec(Psi)) %*% Oneq) / tr(Psi), name="pi")
        ## Equation 7, V: pxm, Lambda: mxm, H: qxm
        B_unstand <- mxAlgebra( V%*%Fpm %*% solve(t(Fpm)%&%Lambda) %*%
                                t(H%*%Fpm), name="B_unstand")
        ## Equation 8
        if (mean.structure) {
            ## beta0 depends on Alpha and Tau; these matrices exist only with mean structure.
            beta0 <- mxAlgebra( Tau - t(B_unstand) %*% V %*% Alpha, name="beta0")
        }
    } else {
        ## Equation 9
        J <- mxAlgebra( H %*% chol(solve(vec2diag(diag2vec(Lambda)))), name="J")
        pi <- mxAlgebra( (t(J*J) %*% Oneq) / tr(Psi), name="pi")
        ## Equation 12, V: pxm, Lambda: mxm, H: qxm
        B_unstand <- mxAlgebra( solve(Dx) %*% V%*%Fpm %*% solve(t(Fpm)%&%Lambda) %*%
                                t(H%*%Fpm) %*% Dy, name="B_unstand")
        ## Equation 13
        if (mean.structure) {
            ## beta0 depends on Alpha and Tau; these matrices exist only with mean structure.
            beta0 <- mxAlgebra( Dy %*% Tau - Dy %*% H%*%Fpm %*% solve(t(Fpm)%&%Lambda) %*%
                                t(V%*%Fpm) %*% V %*% Alpha, name="beta0")
        }
    }

    ## Unbounded value of logit pi
    logitpi <- mxAlgebra( log(pi/(1-pi)), name="logitpi")

    ## Matrix for cumulative sum
    K <- matrix(1, nrow=p, ncol=p)
    K[upper.tri(K)] <- 0
    K <- mxMatrix("Full", nrow=p, ncol=p, free=FALSE, values=K, name="K")

    ## Cumulative sum of pis
    pi_cum <- mxAlgebra(K %*% pi, name="pi_cum")

    if (mean.structure) {
        mx.model <- mxModel(mx.model, J, pi, K, pi_cum, B_unstand, beta0, Fpm, logitpi)
    } else {
        ## Keep beta0 out of the model when mean structure is absent.
        mx.model <- mxModel(mx.model, J, pi, K, pi_cum, B_unstand, Fpm, logitpi)
    }
  
    if (extraTries==0) {
        mx.fit <- suppressMessages(mxRun(mx.model, ...))
    } else {
        mx.fit <- mxTryHard(mx.model, extraTries = extraTries, ...) 
        ## Run it one more time to minimize error
        mx.fit <- suppressMessages(mxRun(mx.fit))
    }

   
    #### Constraints for checking
    ## Diagonals should be either 1 or 0
    Constraint1 <- c(mxEval(vech(V%*%t(V)), mx.fit))

    ## Diagonals should be 1
    Constraint2 <- c(mxEval(diag2vec(V %&% Lambda), mx.fit))

    ## Labels for PCs
    PC_vars <- paste0("PC", seq_len(length(X_vars)))
    PC_cum_vars <- paste0("Sum to PC", seq_len(length(X_vars)))
    
    ## V matrix
    V_est <- mxEval(V, mx.fit)
    V_SE <- suppressMessages(mxSE(V, mx.fit))
    dimnames(V_est) <- dimnames(V_SE) <- list(X_vars, PC_vars)
    
    ## Lambda matrix
    Lambda_est <- mxEval(Lambda, mx.fit)
    Lambda_SE <- suppressMessages(mxSE(Lambda, mx.fit))
    dimnames(Lambda_est) <- dimnames(Lambda_SE) <- list(PC_vars, PC_vars)
    
    ## H matrix
    H_est <- mxEval(H, mx.fit)
    H_SE <- suppressMessages(mxSE(H, mx.fit))
    dimnames(H_est) <- dimnames(H_SE) <- list(Y_vars, PC_vars)

    ## Psi matrix
    Psi_est <- mxEval(Psi, mx.fit)
    Psi_SE <- suppressMessages(mxSE(Psi, mx.fit))
    dimnames(Psi_est) <- dimnames(Psi_SE) <- list(Y_vars, Y_vars)

    ## Unbounded value of logit pi matrix
    ## After calculating the CI of logit pi, it is converted back to the CI of pi
    logitpi_est <- mxEval(logitpi, mx.fit)
    logitpi_SE <- suppressMessages(mxSE(logitpi, mx.fit))
    logitpi_lci <- logitpi_est - 1.96*logitpi_SE
    logitpi_uci <- logitpi_est + 1.96*logitpi_SE

    pi_est <- mxEval(pi, mx.fit)
    pi_lci <- exp(logitpi_lci)/(1+exp(logitpi_lci))
    pi_uci <- exp(logitpi_uci)/(1+exp(logitpi_uci))
    pi <- cbind(pi_est, pi_lci, pi_uci)
    dimnames(pi) <- list(PC_vars, c("Estimate", "Lower 95% CI", "Upper 95% CI"))

    ## Cumulative pi matrix
    pi_cum_est <- mxEval(pi_cum, mx.fit)
    pi_cum_SE <- suppressMessages(mxSE(pi_cum, mx.fit))
    pi_cum_lci <- pi_cum_est - 1.96*pi_cum_SE
    pi_cum_uci <- pi_cum_est + 1.96*pi_cum_SE
    pi_cum <- cbind(pi_cum_est, pi_cum_SE, pi_cum_lci, pi_cum_uci)
    dimnames(pi_cum) <- list(PC_cum_vars, c("Estimate", "SE", "Lower 95% CI", "Upper 95% CI"))
    
    if (pca=="COR") {
        Dx_est <- mxEval(Dx, mx.fit)
        Dx_SE <- suppressMessages(mxSE(Dx, mx.fit))
        dimnames(Dx_est) <- dimnames(Dx_SE) <- list(X_vars, X_vars)

        Dy_est <- mxEval(Dy, mx.fit)
        Dy_SE <- suppressMessages(mxSE(Dy, mx.fit))
        dimnames(Dy_est) <- dimnames(Dy_SE) <- list(Y_vars, Y_vars)    
    } else {
        Dx_est <- Dx_SE <- NULL
        Dy_est <- Dy_SE <- NULL
    }
    
    if (mean.structure) {    
        ## Alpha vector
        Alpha_est <- mxEval(t(Alpha), mx.fit)
        Alpha_SE <- suppressMessages(mxSE(t(Alpha), mx.fit))
        colnames(Alpha_est) <- colnames(Alpha_SE) <- PC_vars

        ## Tau vector
        Tau_est <- mxEval(t(Tau), mx.fit)
        Tau_SE <- suppressMessages(mxSE(t(Tau), mx.fit))
        colnames(Tau_est) <- colnames(Tau_SE) <- Y_vars       
    } else {
        Alpha_est <- Alpha_SE <- NULL
        Tau_est <- Tau_SE <- NULL
    }

    B_unstand_est <- mxEval(B_unstand, mx.fit)
    B_unstand_SE <- suppressMessages(mxSE(B_unstand, mx.fit))
    dimnames(B_unstand_est) <- dimnames(B_unstand_SE) <- list(X_vars, Y_vars)

    if (mean.structure) {
        beta0_est <- mxEval(t(beta0), mx.fit)
        beta0_SE <- suppressMessages(mxSE(t(beta0), mx.fit))
        colnames(beta0_est) <- colnames(beta0_SE) <- Y_vars
    } else {
        beta0_est <- NULL
        beta0_SE <- NULL
    }
           
    out <- list(mx.fit=mx.fit, mx.model=mx.model, Constraint1=Constraint1,
                Constraint2=Constraint2, V_est=V_est, Lambda_est=Lambda_est,
                H_est=H_est, Psi_est=Psi_est, Alpha_est=Alpha_est,
                Tau_est=Tau_est, Dx_est=Dx_est, Dy_est=Dy_est,
                V_SE=V_SE, Lambda_SE=Lambda_SE, H_SE=H_SE, Psi_SE=Psi_SE,
                Alpha_SE=Alpha_SE, Tau_SE=Tau_SE, Dx_SE=Dx_SE, Dy_SE=Dy_SE,
                B_unstand_est=B_unstand_est, B_unstand_SE=B_unstand_SE,
                beta0_est=beta0_est, beta0_SE=beta0_SE, pi=pi, pi_cum=pi_cum,
                pca=pca, pc_select=pc_select)

    class(out) <- "MPCR"  
    out
}


#' Print Method for MPCR Objects
#'
#' Print method for `MPCR` objects.
#'
#' @param x An object returned from the class of `MPCR`.
#' @param digits Number of digits in printing the matrices. The default is 4.
#' @param ... Unused.
#'
#' @return No return value, called for side effects
#'
#' @author Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
#'
#' @method print MPCR
#' @export
print.MPCR <- function(x, digits=4, ...) {
    if (!is.element("MPCR", class(x)))
        stop("\"x\" must be an object of class \"MPCR\".")

    if (x$pca=="COV") {
        cat("\nMPCR: Analysis of covariance matrix.\n")
    } else {
        cat("\nMPCR: Analysis of correlation matrix.\n")
    }
    
    cat("\nPlease check the constraints before interpreting the results.\n")
    cat("Constraint 1. The following values should be either 0 or 1:\n", round(x$Constraint1, 6), ".\n")
    if (x$pca=="COR") cat("Constraint 2. The following values should be 1: ", round(x$Constraint2, 6), ".\n")

    cat("\nThe PC used to construct beta0 and B_unstand are:", x$pc_select, ".\n")
    
    cat("\nV matrix:\n")
    .mprint(x$V_est, digits=digits)
    cat("\nV matrix (SE):\n")
    .mprint(x$V_SE, digits=digits)
  
    cat("\nLambda matrix:\n")
    .mprint(x$Lambda_est, digits=digits)
    cat("\nLambda matrix (SE):\n")
    .mprint(x$Lambda_SE, digits=digits)
  
    cat("\nH matrix:\n")
    .mprint(x$H_est, digits=digits)
    cat("\nH matrix (SE):\n")
    .mprint(x$H_SE, digits=digits)

    cat("\nPsi matrix:\n")
    .mprint(x$Psi_est, digits=digits)
    cat("\nPsi matrix (SE):\n")
    .mprint(x$Psi_SE, digits=digits)
    
    if (!is.null(x$Alpha_est)) {
        cat("\nAlpha vector:\n")
        .mprint(x$Alpha_est, digits=digits)
        cat("\nAlpha vector (SE):\n")
        .mprint(x$Alpha_SE, digits=digits)      
    }

    if (!is.null(x$Dx_est)) {
        cat("\nDx matrix:\n")
        .mprint(x$Dx_est, digits=digits)
        cat("\nDx matrix (SE):\n")
        .mprint(x$Dx_SE, digits=digits)      
    }

    if (!is.null(x$Dy_est)) {
        cat("\nDy matrix:\n")
        .mprint(x$Dy_est, digits=digits)
        cat("\nDy matrix (SE):\n")
        .mprint(x$Dy_SE, digits=digits)      
    }
        
    if (!is.null(x$Tau_est)) {
        cat("\nTau vector:\n")
        .mprint(x$Tau_est, digits=digits)
        cat("\nTau vector (SE):\n")
        .mprint(x$Tau_SE, digits=digits)      
    }


    if (!is.null(x$beta0_est)) {
        cat("\nbeta0 vector:\n")
        .mprint(x$beta0_est, digits=digits)
        cat("\nbeta0 vector (SE):\n")
        .mprint(x$beta0_SE, digits=digits)
    }
    
    cat("\nB_unstand matrix:\n")
    .mprint(x$B_unstand_est, digits=digits)
    cat("\nB_unstand matrix (SE):\n")
    .mprint(x$B_unstand_SE, digits=digits)

    cat("\npi vector:\n")
    .mprint(x$pi, digits=digits)

    cat("\nCumulative pi vector:\n")
    .mprint(x$pi_cum, digits=digits)     
}
