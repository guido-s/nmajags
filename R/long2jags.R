#' Transform data from long arm-based to JAGS
#' format
#' 
#' @description
#' This function transforms data from long arm-based format, i.e., two
#' rows for a pairwise comparison, to a list format for Bayesian
#' network meta-analysis using JAGS.
#' 
#' @param studlab A vector with study labels or an object created with
#'   \code{\link[meta]{longarm}}.
#' @param treat A vector with treatment labels (mandatory if
#'   \code{studlab} is a vector).
#' @param event Number of events.
#' @param n Number of observations.
#' @param mean Estimated mean.
#' @param sd Standard deviation.
#' @param data An optional data frame containing the study
#'   information.
#' @param reference.group Reference treatment.
#' @param covar Covariable to use in network meta-regression.
#' @param func A character string specifying function to summarize
#'   data for additional covariate. Either, "min", "max", "mean", or
#'   ""; can be abbreviated.
#' 
#' @details
#' This function transforms data from long arm-based format, i.e., two
#' rows for a pairwise comparison, to a list format for Bayesian
#' network meta-analysis using JAGS.
#' 
#' At the moment, the function can be used to transform data with a
#' binary or continuous outcome. The following arguments are mandatory
#' if argument \code{studlab} is not an R object created with
#' \code{\link[meta]{longarm}}. Arguments \code{studlab} and
#' \code{treat} must be provided to identify studies and treatments
#' and, depending on the outcome, the following additional arguments:
#' 
#' \itemize{
#' \item \code{event}, \code{n} (binary outcome);
#' \item \code{n}, \code{mean}, \code{sd} (continuous outcome).
#' }
#' 
#' @return
#' A list of class \code{jagsdata}.
#' 
#' @author Georgia Salanti \email{georgia.salanti@@ispm.unibe.ch},
#'   Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netjags}}, \code{\link[meta]{longarm}},
#'   \code{\link[netmeta]{netmeta}}
#' 
#' @keywords datagen
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
#' 
#' # Transform long-arm based data into a list suitable for JAGS
#' # analysis, with placebo as reference treatment
#' #
#' dat.jags <- long2jags(study, trt,
#'                       event = stroke, n = total,
#'                       data = Dogliotti2014, 
#'                       reference = "pl")
#' 
#' # Run JAGS and create JAGS object
#' #
#' n1.jags <- netjags(dat.jags)
#' 
#' # Extract lnOR from BUGS output
#' #
#' trts <- attr(n1.jags, "trts")
#' TE.random <- n1.jags$BUGSoutput$mean$LOR.full
#' rownames(TE.random) <- trts
#' colnames(TE.random) <- trts
#' 
#' # Extract se(lnOR) from BUGS output
#' #
#' seTE.random <- n1.jags$BUGSoutput$sd$LOR.full
#' rownames(seTE.random) <- trts
#' colnames(seTE.random) <- trts
#' 
#' # Conduct frequentist NMA with tau from Bayesian NMA
#' #
#' tau.jags <- as.numeric(n1.jags$BUGSoutput$mean$tau)
#' n1.iv.jags <-
#'   suppressWarnings(netmeta(p1, tau.preset = tau.jags))
#' 
#' # Print treatment matrix for frequentist NMA
#' #
#' round(n1.iv$TE.random, 3)[trts, trts]
#' 
#' # Print treatment matrix for Bayesian NMA
#' #
#' round(TE.random, 3)
#' 
#' # Print treatment matrix for frequentist NMA
#' # (using tau from Bayesian NMA)
#' #
#' round(n1.iv.jags$TE.random, 3)[trts, trts]
#' 
#' # Print standard errors from frequentist NMA
#' #
#' round(n1.iv$seTE.random, 3)[trts, trts]
#' 
#' # Print standard errors from Bayesian NMA
#' #
#' round(seTE.random, 3)
#' 
#' # Print standard errors from frequentist NMA
#' # (using tau from Bayesian NMA)
#' #
#' round(n1.iv.jags$seTE.random, 3)[trts, trts]
#' 
#' # Print square root of between-study heterogeneity tau2
#' #
#' tau.jags
#' n1.iv$tau
#' }
#' 
#' @export long2jags


long2jags <- function(studlab, treat,
                      event, n, mean, sd,
                      data = NULL, reference.group = 1,
                      covar, func = "max") {
  
  
  ##
  ##
  ## (1) Read data
  ##
  ##
  nulldata <- is.null(data)
  if (nulldata)
    data <- sys.frame(sys.parent())
  mf <- match.call()
  ##
  if (missing(studlab))
    stop("Mandatory argument 'studlab' missing.")
  else
    studlab <- eval(mf[[match("studlab", names(mf))]],
                    data, enclos = sys.frame(sys.parent()))
  ##
  ## R object created with longarm()
  ##
  if (is.data.frame(studlab) & !is.null(attr(studlab, "longarm"))) {
    type <- attr(studlab, "type")
    ##
    treat <- studlab$treat
    n <- studlab$n
    ##
    if (type == "binary") {
      event <- studlab$events
      mean <- NULL
      sd <- NULL
    }
    else if (type == "continuous") {
      event <- NULL
      mean <- studlab$mean
      sd <- studlab$sd
    }
    else
      stop("Function cannot be used with count or generic outcome.")
    ##
    studlab <- studlab$studlab
  }
  else {
    if (missing(treat))
      stop("Argument 'treat' mandatory.")
    else
      treat <- eval(mf[[match("treat", names(mf))]],
                    data, enclos = sys.frame(sys.parent()))
    ##
    ## Catch event, n, mean, sd, covar from data:
    ##
    event <- eval(mf[[match("event", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
    n <- eval(mf[[match("n", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
    mean <- eval(mf[[match("mean", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
    sd <- eval(mf[[match("sd", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  }
  ##
  missing.covar <- missing(covar)
  if (!missing.covar)
    covar <- eval(mf[[match("covar", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
  ##
  if (!is.null(event))
    chknumeric(event)
  ##
  if (!is.null(n))
    chknumeric(n)
  ##
  if (!is.null(mean))
    chknumeric(mean)
  ##
  if (!is.null(sd))
    chknumeric(sd)
  ##
  if (!is.null(event) & !is.null(n) &
      is.null(mean) & is.null(sd))
    type <- "binary"
  else if (!is.null(n) & !is.null(mean) & !is.null(sd))
    type <- "continuous"
  else
    stop("Type of outcome unclear. Please provide the necessary ",
         "information:\n  ",
         "- event, n (binary outcome)\n  ",
         "- n, mean, sd (continuous outcome).")
  
  
  ##
  ##
  ## (2) Check length of variables
  ##
  ##
  Nobs <- length(studlab)
  ##
  chklength(treat, Nobs, name = "studlab")
  ##
  if (!is.null(event))
    chklength(event, Nobs, name = "studlab")
  ##
  if (!is.null(n))
    chklength(n, Nobs, name = "studlab")
  ##
  if (!is.null(mean))
    chklength(mean, Nobs, name = "studlab")
  ##
  if (!is.null(sd))
    chklength(sd, Nobs, name = "studlab")
  ##
  if (!missing.covar)
    chklength(covar, Nobs, name = "studlab")
  ##
  chklength(reference.group, 1)
  
  
  ##
  ##
  ## (3) Data set
  ##
  ##
  if (nulldata) {
    data <- data.frame(.studlab = studlab, .treat = treat, .n = n)
    ##
    if (type == "binary")
      data$.event <- event
    else if (type == "continuous") {
      data$.mean <- mean
      data$.sd <- sd
    }
    ##
    if (!missing.covar)
      data$.covar <- covar
  }
  else {
    data$.studlab <- studlab
    data$.treat <- treat
    data$.event <- event
    data$.n <- n
    data$.mean <- mean
    data$.sd <- sd
    if (!missing.covar)
      data$.covar <- covar
  }
  
  
  ##
  ##
  ## (4) Prepare data for JAGS list
  ##
  ##
  trts <- sort(unique(treat))
  id.trts <- seq_along(trts)
  ##
  reference.group <- setref(reference.group, trts)
  if (trts[1] != reference.group)
    trts <- c(reference.group, trts[!(trts == reference.group)])
  ##
  id.studlab <- as.numeric(as.factor(studlab))
  id.treat <- as.numeric(factor(treat, trts, id.trts))
  ##
  o <- order(id.studlab, id.treat)
  ##
  id.studlab <- id.studlab[o]
  id.treat <- id.treat[o]
  ##
  studlab <- studlab[o]
  treat <- treat[o]
  n <- n[o]
  ##
  if (type == "continuous") {
    mean <- mean[o]
    sd <- sd[o]
  }
  else if (type == "binary")
    event <- event[o]
  ##
  studs <- unique(studlab)
  
  
  ##
  ##
  ## (5) Generate JAGS list
  ##
  ##
  extract <- function(x, func) {
    setchar(func, c("max", "min", "mean"))
    ##
    if (func == "max")
      apply(x, 1, max, na.rm = TRUE)
    else if (func == "min")
      apply(x, 1, min, na.rm = TRUE)
    else if (func == "mean")
      apply(x, 1, mean, na.rm = TRUE)
    else
      stop("Admissible values for argument 'func': ",
           "\"max\", \"min\", \"mean\" or \"\".")
  }
  ##
  asMat <- function(x, id1, id2, rownames, colnames) {
    res <- matrix(NA, nrow = length(rownames), ncol = length(colnames))
    ##
    for (i in seq_along(x)) {
      res[id1[i], id2[i]] <- x[i]
    }
    ##
    rownames(res) <- if (is.numeric(rownames))
                       as.character(rownames)
                     else
                       rownames
    colnames(res) <- if (is.numeric(colnames))
                       as.character(colnames)
                     else
                       colnames
    ##
    res
  }
  ##
  n.k <- table(id.studlab)
  T <- asMat(id.treat,
             id.studlab, unlist(sapply(n.k, seq)),
             studs, seq_len(max(table(id.studlab))))
  T[is.na(T)] <- 666
  T <- t(apply(T, 1, sort))
  T[T == 666] <- NA
  ##
  res <- list(k = length(studs), n = length(trts), n.k = n.k, T = T)
  ##
  if (type == "continuous") {
    res[["Y"]] <- asMat(mean, id.studlab, id.treat, studs, trts)
    ##
    res[["Pr"]] <-
      asMat(1 / (sd / sqrt(n))^2, id.studlab, id.treat, studs, trts)
    ##
    res[["pooled.sd"]] <-
      sqrt(tapply(n * sd^2, id.studlab, sum)) /
      sqrt(tapply(n, id.studlab, sum) - n.k)
  }
  else if (type == "binary") {
    res[["E"]] <- asMat(event, id.studlab, id.treat, studs, trts)
    ##
    res[["N"]] <- asMat(n, id.studlab, id.treat, studs, trts)
  }
  ##
  if (!missing(covar)) {
    res[["func"]] <- func
    res[["Covar"]] <- asMat(covar, id.studlab, id.treat, studs, trts)
    ##
    if (func != "")
      res[["Covar"]] <- extract(res[["Covar"]], func)
  }
  ##
  res[["trts"]] <- trts
  res[["data"]] <- data
  ##  
  class(res) <- "jagsdata"
  attr(res, "type") <- type
  ##
  res
}
