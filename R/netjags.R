#' Bayesian network meta-analysis using JAGS
#' 
#' @description
#' Network meta-analysis for objects of class \code{jagsdata}. This is
#' a wrapper function for the R function \code{\link[R2jags]{jags}} in
#' the R package \bold{R2jags} (Su & Yajima 2021).
#' 
#' @param x An object of class \code{jagsdata}.
#' @param inits ...
#' @param parameters.to.save ...
#' @param n.chains ...
#' @param n.iter ...
#' @param n.burnin ...
#' @param DIC ...
#' @param model.file ...
#' @param quiet A logical indicating whether information on the
#'   estimation progress should be printed.
#' @param \dots Additional arguments (passed on to
#'   \code{\link[R2jags]{jags}})
#' 
#' @details
#' Network meta-analysis for objects of class \code{jagsdata} which
#' can be created with \code{\link{long2jags}}. This is a wrapper
#' function for the R function \code{\link[R2jags]{jags}} in the R
#' package \bold{R2jags} (Su & Yajima 2021).
#' 
#' @author Georgia Salanti \email{georgia.salanti@@ispm.unibe.ch},
#'   Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{long2jags}}, \code{\link[netmeta]{netmeta}}
#' 
#' @references
#' Yu-Sung Su and Masanao Yajima, (2021):
#' R2jags: Using R to Run 'JAGS'.
#' https://CRAN.R-project.org/package=R2jags
#' 
#' @examples
#' \dontrun{
#' library(netmeta)
#' data(Dogliotti2014)
#' Dogliotti2014$trt <- treats(Dogliotti2014$treatment, nchar = 4)
#' #
#' p1 <- pairwise(trt, stroke, total, studlab = study,
#'                data = Dogliotti2014, sm = "OR")
#' 
#' # Standard frequentist random effects NMA model
#' #
#' n1.iv <- netmeta(p1)
#' 
#' # Standard Bayesian random effects NMA model
#' #
#' library(R2jags)
#' #
#' # Transform long-arm based data into a list suitable for JAGS
#' # analysis, with first treatment (here: Apixaban) as reference
#' # treatment, i.e., t = 1
#' #
#' dat.jags <- long2jags(id, trt,
#'                       event = stroke, n = total,
#'                       data = Dogliotti2014, 
#'                       reference = n1.iv$trts[1])
#' dat.jags
#' #
#' # Run JAGS and create JAGS object (run 10.000 iterations and 1.000
#' # burn-in)
#' #
#' n1.jags <- netjags(dat.jags)
#' #
#' # Print results
#' #
#' print(n1.jags, digits = 2)
#' 
#' # Extract lnOR from Bayesian NMA
#' #
#' TE.random <- n1.jags$BUGSoutput$mean$LOR.full
#' rownames(TE.random) <- rownames(n1.iv$TE.random)
#' colnames(TE.random) <- colnames(n1.iv$TE.random)
#' 
#' # Extract se(lnOR) from Bayesian NMA
#' #
#' seTE.random <- n1.jags$BUGSoutput$sd$LOR.full
#' rownames(seTE.random) <- rownames(n1.iv$seTE.random)
#' colnames(seTE.random) <- colnames(n1.iv$seTE.random)
#' 
#' # Conduct frequentist NMA with tau from Bayesian NMA
#' #
#' tau.jags <- as.numeric(n1.jags$BUGSoutput$mean$tau)
#' n1.iv.jags <-
#'   suppressWarnings(netmeta(p1, tau.preset = tau.jags))
#' 
#' # Print treatment matrix for frequentist NMA
#' #
#' round(n1.iv$TE.random, 3)
#' 
#' # Print treatment matrix for Bayesian NMA
#' #
#' round(TE.random, 3)
#' 
#' # Print treatment matrix for frequentist NMA
#' # (using tau from Bayesian NMA)
#' #
#' round(n1.iv.jags$TE.random, 3)
#' 
#' # Print standard errors from frequentist NMA
#' #
#' round(n1.iv$seTE.random, 3)
#' 
#' # Print standard errors from Bayesian NMA
#' #
#' round(seTE.random, 3)
#' 
#' # Print standard errors from frequentist NMA
#' # (using tau from Bayesian NMA)
#' #
#' round(n1.iv.jags$seTE.random, 3)
#' 
#' # Print square root of between-study heterogeneity tau2
#' #
#' tau.jags
#' n1.iv$tau
#' }
#' 
#' @export netjags


netjags <- function(x, inits = NULL,
                    parameters.to.save = NULL,
                    n.chains = 2, n.iter = 10000,
                    n.burnin = 1000, DIC = FALSE,
                    model.file = NULL,
                    quiet = FALSE,
                    ...) {
  
  chkclass(x, "jagsdata")
  type <- attr(x, "type")
  ##
  chknumeric(n.chains, min = 1, length = 1)
  chknumeric(n.iter, min = 1, length = 1)
  chknumeric(n.burnin, min = 1, length = 1)
  chklogical(DIC)
  chklogical(quiet)
  ##
  if (is.null(model.file))
    if (type == "binary")
      model.file <-
        system.file("model", "nma.bin.or.txt", package = "nmajags")
    else if (type == "continuous")
      model.file <-
        system.file("model", "nma.cont.smd.txt", package = "nmajags")
    else
      stop("Argument 'model.file' must be provided.")
  ##
  if (is.null(parameters.to.save))
    if (type == "binary")
      parameters.to.save <- c("LOR.full", "tau")
    else if (type == "continuous")
      parameters.to.save <- c("SMD.full", "tau")
    else
      stop("Argument 'parameters.to.save' must be provided.")
  ##
  x$data <- NULL
  ##
  res <- jags(data = x, inits = inits,
              parameters.to.save = parameters.to.save,
              n.chains = n.chains, n.iter = n.iter,
              n.burnin = n.burnin, DIC = DIC,
              model.file = model.file, quiet = quiet,
              ...)
  ##
  res
}
