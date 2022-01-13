#' Transform data from long arm-based to JAGS
#' format
#' 
#' @description
#' This function transforms data from long arm-based format, i.e., two
#' rows for a pairwise comparison, to a list format for Bayesian
#' network meta-analysis using JAGS.
#' 
#' @param studlab A vector with study labels (mandatory).
#' @param treat A vector with treatment labels (mandatory).
#' @param event Number of events.
#' @param n Number of observations.
#' @param mean Estimated mean.
#' @param sd Standard deviation.
#' @param data An optional data frame containing the study
#'   information.
#' @param reference.group Reference treatment.
#' @param addvar Additional covariable to use in network
#'   meta-regression.
#' @param func A character string specifying function to summarize
#'   data for additional covariate. Either, "min", "max" or "mean",
#'   can be abbreviated.
#' @param quiet A logical indicating whether information on changes in
#'   treatment labels should be printed.
#' 
#' @details
#' This function transforms data from long arm-based format, i.e., two
#' rows for a pairwise comparison, to a list format for Bayesian
#' network meta-analysis using JAGS.
#' 
#' The function can be used to transform data with a binary or
#' continuous outcome. Arguments \code{studlab} and \code{treat} are
#' mandatory to identify studies and treatments and, depending on the
#' outcome, the following additional arguments are mandatory:
#' 
#' \itemize{
#' \item \code{event}, \code{n} (binary outcome);
#' \item \code{n}, \code{mean}, \code{sd} (continuous outcome).
#' }
#' 
#' @return
#' A list of class \code{jagsdata} readable in JAGS.
#' 
#' @author Georgia Salanti \email{georgia.salanti@@ispm.unibe.ch},
#'   Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netjags}}, \code{\link[netmeta]{netmeta}}
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
#' # Run JAGS and create JAGS object (run 10.000 iterations and
#' # 1.000 burn-in) monitor the nodes LOR and tau
#' #
#' n1.jags <- netjags(dat.jags)
#' #
#' # Print results
#' #
#' print(n1.jags, digits = 2)
#' 
#' # Extract lnOR from Bayesian NMA
#' #
#' TE.random <- n1.jags$BUGSoutput$mean$lnOR.full
#' rownames(TE.random) <- rownames(n1.iv$TE.random)
#' colnames(TE.random) <- colnames(n1.iv$TE.random)
#' 
#' # Extract selnOR from Bayesian NMA
#' #
#' seTE.random <- n1.jags$BUGSoutput$sd$lnOR.full
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
#' @export long2jags


long2jags <- function(studlab, treat,
                      event, n, mean, sd,
                      data = NULL, reference.group = 1,
                      addvar, func = "max",
                      quiet = FALSE) {
  
  
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
    stop("Argument 'studlab' mandatory.")
  if (missing(treat))
    stop("Argument 'treat' mandatory.")
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  treat <- eval(mf[[match("treat", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  ##
  ##
  ## Catch event, n, mean, sd, addvar from data:
  ##
  event <- eval(mf[[match("event", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  n <- eval(mf[[match("n", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  mean <- eval(mf[[match("mean", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  sd <- eval(mf[[match("sd", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  missing.addvar <- missing(addvar)
  if (!missing.addvar)
    addvar <- eval(mf[[match("addvar", names(mf))]],
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
  chklogical(quiet)
  
  
  ##
  ##
  ## (2) Check length of variables
  ##
  ##
  k.all <- length(studlab)
  chklength(treat, k.all, name = "studlab")
  ##
  if (!is.null(event))
    chklength(event, k.all, name = "studlab")
  ##
  if (!is.null(n))
    chklength(n, k.all, name = "studlab")
  ##
  if (!is.null(mean))
    chklength(mean, k.all, name = "studlab")
  ##
  if (!is.null(sd))
    chklength(sd, k.all, name = "studlab")
  ##
  if (!missing.addvar)
    chklength(addvar, k.all, name = "studlab")
  ##
  chklength(reference.group, 1)
  
  
  ##
  ##
  ## (3) Data set
  ##
  ##
  if (nulldata) {
    data <-
      data.frame(.studlab = studlab,
                 .treat = treat,
                 .event = event, .n = n,
                 .mean = mean, .sd = sd)
    if (!missing.addvar)
      data$.addvar <- addvar
  }
  else {
    data$.studlab <- studlab
    data$.treat <- treat
    data$.event <- event
    data$.n <- n
    data$.mean <- mean
    data$.sd <- sd
    if (!missing.addvar)
      data$.addvar <- addvar
  }
  
  
  ##
  ##
  ## (4) Prepare data for JAGS list
  ##
  ##
  s.id <- as.numeric(as.factor(studlab))
  t.id <- as.numeric(as.factor(treat))
  ##
  ns <- length(unique(s.id))
  nt <- length(unique(treat))
  ##
  na <- table(s.id)
  s.id <- rep(1:ns, table(s.id))
  ##  
  if (!identical(t.id, treat)) {
    cat("Note: the treatments have been renamed as follows:\n")
    out <- cbind.data.frame(`old names` = sort(unique(treat)),
                            `new names` = sort(unique(t.id)))
    reference.group <-
      sort(unique(t.id))[sort(unique(treat)) == reference.group]
    ##
    if (!quiet) {
      out <- as.matrix(out)
      rownames(out) <- rep("", nrow(out))
      prmatrix(out, quote = FALSE, right = TRUE)
    }
  }
  ##  
  maxnrofarms <- max(table(s.id))
  n.arms <- length(s.id)
  armsenumerate <- unlist(sapply(na, seq))
  ##
  tmat <- matrix(666.666, nrow = ns, ncol = maxnrofarms)
  for (i in 1:n.arms)
    tmat[s.id[i], armsenumerate[i]] <- t.id[i]
  ##
  tmat <- t(apply(tmat, 1, sort))
  tmat[tmat == 666.666] <- NA
  ##
  nmat <- matrix(-99, nrow = ns, ncol = nt)
  for (i in 1:n.arms) {
    nmat[s.id[i], t.id[i]] <- n[i]
  }
  nmat[nmat == -99] <- NA
  ##
  if (!missing(addvar)) {
    vmat <- matrix(NA, nrow = ns, ncol = nt)
    for (i in 1:n.arms) {
      vmat[s.id[i], t.id[i]] <- addvar[i]
    }
    vmat[is.na(vmat)] <- NA
  }
  
  
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
           "\"max\", \"min\", \"mean\".")
  }
  ##
  if (type == "continuous") {
    ymat <- pmat <- matrix(-666, nrow = ns, ncol = nt)
    ##
    prec <- 1 / (sd / sqrt(n))^2
    ##
    for (i in 1:n.arms) {
      ymat[s.id[i], t.id[i]] <- mean[i]
      pmat[s.id[i], t.id[i]] <- prec[i]
    }
    ##
    ymat[ymat == -666] <- NA
    pmat[pmat == -666] <- NA
    ##
    ## Calculate pooled SD for SMD
    ##
    nominator <- sqrt(tapply(n * sd^2, s.id, sum))
    denominator <- sqrt(tapply(n, s.id, sum) - na)
    pooled.sd <- nominator / denominator
    ##
    res <- list(ns = ns, nt = nt, na = na, t = tmat,
                y = ymat, pooled.sd = pooled.sd,
                prec = pmat, ref = reference.group)
  }
  else if (type == "binary") {
    emat <- matrix(-666, nrow = ns, ncol = nt)
    ##
    for (i in 1:n.arms) {
      emat[s.id[i], t.id[i]] <- event[i]
    }
    emat[emat == -666] <- NA
    ##
    res <- list(ns = ns, nt = nt, na = na, t = tmat,
                r = emat, n = nmat, ref = reference.group)
  }
  
  
  if (!missing(addvar))
    res[["addvar"]] <- extract(vmat, func)
  ##
  res[["data"]] <- data
  ##  
  class(res) <- "jagsdata"
  attr(res, "type") <- type
  ##
  res
}
