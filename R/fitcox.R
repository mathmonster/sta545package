#' Cox proportional hazards model
#'
#' \code{fitcox} fits a Cox proportional hazards model.
#'
#' @param formula A formula object.
#' @param delta Censoring indicator.
#' @param data A data frame.
#' @return An object of class \code{fitcox}.
#' @details The object returned has the attributes:
#' \describe{
#'    \item{formula}{The formula used to generate the model}
#'    \item{delta}{The column containing the censoring indicator}
#'    \item{n}{The total number of observations}
#'    \item{n.events}{The number of events i.e. uncensored observations}
#'    \item{n.times}{The number of distinct event times}
#'    \item{event.times}{The distinct event times}
#'    \item{event.counts}{The number of events taking place at each event time}
#'    \item{h0}{The baseline hazard at each event time}
#' }
#' @examples
#' fitcox(remission ~ sample, delta="censor", data=leukemia)
#' @seealso \code{\link{summary.fitcox}}
#' @importFrom stats model.matrix optim pnorm qnorm pchisq
#' @export

fitcox <- function(formula, delta, data) {
  n <- nrow(data)
  response.var <- all.vars(formula)[1]
  # Extract the variables needed to fit the formula
  responses <- data[[response.var]]
  deltas <- data[[delta]]
  # Observations for which delta = 1 i.e. events

  # Extract the event times
  all.times <- data[[response.var]]
  events <- all.times[deltas==1]
  event.times <- sort(unique(events))
  k <- length(event.times)

  # Event sets are the sets of non-censored observations corresponding to the
  # event times
  event.sets <- lapply(event.times,
                       function(t) model.matrix(formula,
                                                data=data[responses == t &
                                                            deltas == 1, ]))

  event.counts <- vapply(event.sets, nrow, numeric(1))

  # Number of coefficients is determined by the factor levels found by
  # model.matrix
  p <- ncol(event.sets[[1]]) - 1

  # Risk set differences are the sets of all observations with
  # t(i) <= t < t(i + 1). These are the observations at risk at each event time,
  # but not at the next
  risk.set.diffs <- list()
  for (i in 1:(k-1)) {
    t1 <- event.times[i]
    t2 <- event.times[i + 1]
    risk.set.diffs[[i]] <- model.matrix(formula, data=data[responses >= t1 &
                                                             responses < t2, ])
  }
  # The last risk set difference contains the risk set for the final event time
  risk.set.diffs[[k]] <- model.matrix(formula,
                                      data=data[responses >= event.times[k], ])

  # Use optim to find the MLEs for the regression coefficients
  beta <- optim(as.matrix(rep(0, p)), partial.log.likelihood, gr=partial.score,
                event.sets, risk.set.diffs,
                method="BFGS", control=list(fnscale=-1))$par

  # Calculate standard errors using the observed information matrix
  information <- partial.information(beta, event.sets, risk.set.diffs)
  beta.se <- as.vector(1/sqrt(diag(information)))
  alpha <- 0.05
  z.lower <- qnorm(alpha/2)
  z.upper <- qnorm(1 - alpha/2)
  beta.lower <- beta + z.lower*beta.se
  beta.upper <- beta + z.upper*beta.se

  # Wald test used for coefficients
  beta.z <- beta/beta.se
  beta.p <- 2*pnorm(-abs(beta.z))

  # Store baseline hazard for making predictions
  baseline <- baseline.hazard(beta, event.counts, risk.set.diffs)

  # Likelihood ratio test
  lrt <- 2*(partial.log.likelihood(beta, event.sets, risk.set.diffs)
            - partial.log.likelihood(rep(0, p), event.sets, risk.set.diffs))
  lrt.p <- 1 - pchisq(lrt, df=1)

  results <- list(formula=formula, delta=delta,
                  n=n, n.events=length(events), n.times=k,
                  parameters=colnames(event.sets[[1]])[-1],
                  event.times=event.times, event.counts=event.counts,
                  h0=baseline,
                  lrt=lrt, lrt.p=lrt.p,
                  beta=beta, beta.se=beta.se,
                  beta.lower=beta.lower, beta.upper=beta.upper,
                  beta.z=beta.z, beta.p=beta.p)
  class(results) <- "fitcox"
  return(results)
}

#' Log partial likelihood for Cox model using Efron's method for ties
#'
#' @param beta A vector of length p, where p is the number of model parameters
#' @param event.sets A list of model matrices corresponding to the event times
#' @param risk.set.diffs A list of risk set differences for the event times
#' @return The log partial likelihood, a scalar

partial.log.likelihood <- function(beta, event.sets, risk.set.diffs) {
  k <- length(event.sets)
  l <- 0
  beta <- c(0, beta)

  # R0 is the running total of exp(z*beta) over the risk set
  R0 <- 0

  # Efron's method is performed over the reversed event times
  # The risk set is the union of all previous risk set differences
  for (i in k:1) {
    event.set <- event.sets[[i]]
    risk.set.diff <- risk.set.diffs[[i]]
    d <- nrow(event.set)
    # R0 is the sum of exp(z*beta) over the risk set
    R0 <- R0 + sum(exp(risk.set.diff%*%beta))
    # E0 is the sum of exp(z*beta) over the event set
    events.beta <- event.set%*%beta
    E0 <- sum(exp(events.beta))
    l <- l + sum(events.beta) - sum(log(R0 - E0*(0:(d - 1))/d))
  }
  return(l)
}

#' Score function for Cox model using Efron's method for ties
#'
#' @param beta A vector of length p, where p is the number of model parameters
#' @param event.sets A list of model matrices corresponding to the event times
#' @param risk.set.diffs A list of risk set differences for the event times
#' @return A vector of partial derivatives for each parameter in the model

partial.score <- function(beta, event.sets, risk.set.diffs) {
  k <- length(event.sets)
  p <- length(beta)

  # First term in u and beta is the intercept: this gets discarded at the end
  u <- rep(0, p + 1)
  beta <- c(0, beta)

  # Running totals of exp(z*beta) and z*exp(z*beta) over the risk sets
  R0 <- 0
  R1 <- rep(0, p + 1)

  # Score function for Efron's method, performed over the reversed event times
  for (i in k:1) {
    event.set <- event.sets[[i]]
    risk.set.diff <- risk.set.diffs[[i]]
    d <- nrow(event.set)

    events.beta <- event.set%*%beta
    risk.set.beta <- risk.set.diff%*%beta
    exp.events.beta <- exp(events.beta)
    exp.risk.set.beta <- exp(risk.set.beta)

    R0 <- R0 + sum(exp.risk.set.beta)
    E0 <- sum(exp.events.beta)
    R1 <- R1 + as.vector(t(risk.set.diff)%*%exp.risk.set.beta)
    E1 <- as.vector(t(event.set)%*%exp.events.beta)
    # First term is the sum of the columns in the event set
    # Second term is the sum over the event set of the log Efron denominator
    u <- u + apply(event.set, 2, sum) -
      Reduce('+', lapply(0:(d-1), function(m) (R1 - m*E1/d)/(R0 - m*E0/d)))
  }
  return(u[-1])
}

#' Information matrix for Cox model using Efron's method for ties
#'
#' @param beta A vector of length p, where p is the number of model parameters
#' @param event.sets A list of model matrices corresponding to the event times
#' @param risk.set.diffs A list of risk set differences for the event times
#' @return A matrix of negative 2nd partial derivatives for each parameter pair

partial.information <- function(beta, event.sets, risk.set.diffs) {
  k <- length(event.sets)
  p <- length(beta)

  # First term in u and beta is the intercept: this gets discarded at the end
  information <- matrix(0, nrow = p + 1, ncol = p + 1)
  beta <- c(0, beta)

  # Totals of exp(z*beta), z*exp(z*beta), and z*z*exp(z*beta) over the risk sets
  R0 <- 0
  R1 <- rep(0, p + 1)
  R2 <- matrix(0, nrow = p + 1, ncol = p + 1)

  # Terms to be summed
  efron.denom <- function(m, d, R0, E0, R1, E1, R2, E2) {
    g0 <- R0 - m*E0/d
    g1 <- R1 - m*E1/d
    g2 <- R2 - m*E2/d
    return((g0*g2 - outer(g1, g1))/(g0^2))
  }

  # Hessian for Efron's method, performed over the reversed event times
  # Notation Ri and Ei denotes sum of ith derivatives for risk and event sets
  for (i in k:1) {
    event.set <- event.sets[[i]]
    risk.set.diff <- risk.set.diffs[[i]]
    d <- nrow(event.set)

    events.beta <- event.set%*%beta
    risk.set.beta <- risk.set.diff%*%beta
    exp.events.beta <- exp(events.beta)
    exp.risk.set.beta <- exp(risk.set.beta)

    R0 <- R0 + sum(exp.risk.set.beta)
    E0 <- sum(exp.events.beta)
    R1 <- R1 + as.vector(t(risk.set.diff)%*%exp.risk.set.beta)
    E1 <- as.vector(t(event.set)%*%exp.events.beta)
    R2 <- R2 + matrix(rowSums(apply(risk.set.diff, 1,
                                    function(z) outer(z, z)*exp(sum(z*beta)))),
                      nrow = p + 1)
    E2 <- matrix(rowSums(apply(event.set, 1,
                               function(z) outer(z, z)*exp(sum(z*beta)))),
                 nrow = p + 1)
    # Sum over m = 0, ..., d - 1
    information <- information + Reduce('+', lapply(0:(d-1), efron.denom,
                                                    d, R0, E0, R1, E1, R2, E2))
  }
  return(matrix(information[-1, -1], nrow=p))
}

#' Baseline hazard for a Cox proportional hazards model
#'
#' @param beta A vector of length p, where p is the number of model parameters
#' @param event.counts A vector of event counts for the unique event times
#' @param risk.set.diffs A list of risk set differences for the event times
#' @return A vector of partial derivatives for each parameter in the model

baseline.hazard <- function(beta, event.counts, risk.set.diffs) {
  k <- length(risk.set.diffs)
  hazard <- rep(NA, k)

  # First term in beta is the intercept: set to zero
  beta <- c(0, beta)

  # Running total of exp(z*beta) over the risk sets
  R.sum <- 0

  # Baseline hazard performed over the reversed event times
  for (i in k:1) {
    risk.set.diff <- risk.set.diffs[[i]]
    R.sum <- R.sum + sum(exp(risk.set.diff%*%beta))
    hazard[i] <- event.counts[i]/R.sum
  }

  return(hazard)
}

#' Extract the coefficients for a fitted Cox model
#'
#' @param model A \code{fitcox} object returned by the \code{fitcox} function.
#' @return A dataframe with one row per coefficient
#' @details Coefficient details include:
#' \describe{
#'    \item{beta}{The value of the coefficient}
#'    \item{beta.se}{Standard error for the coefficient}
#'    \item{beta.lower}{Lower bound of 95\% confidence interval}
#'    \item{beta.upper}{Upper bound of 95\% confidence interval}
#'    \item{z}{Wald statistic i.e. beta/beta.se}
#'    \item{p}{p-value for the Wald statistic}
#' }
#' @examples
#' coef(fitcox(remission ~ sample, delta="censor", data=leukemia))
#' @importFrom stats coef
#' @export

coef.fitcox <- function(model) {
  return(data.frame(beta=model$beta,
                    beta.se=model$beta.se,
                    lower=model$beta.lower,
                    upper=model$beta.upper,
                    z=model$beta.z,
                    p=model$beta.p,
                    row.names=model$parameters))
}

#' Print details of a fitted Cox model
#'
#' @param model A \code{fitcox} object returned by the \code{fitcox} function.
#' @return The object of class \code{fitcox}
#' @examples
#' print(fitcox(remission ~ sample, delta="censor", data=leukemia))
#' @export

print.fitcox <- function(model) {
  print(coef(model), digits=4)
  cat("\n")
  cat("n = ", model$n, ", number of events = ", model$n.events,
  ", distinct event times = ", model$n.times, "\n", sep="")
  cat("Likelihood ratio test: ", model$lrt, ", p = ", model$lrt.p, sep="")
  invisible(model)
}

#' Generate a summary for a fitted Cox model
#'
#' @param model A \code{fitcox} object returned by the \code{fitcox} function.
#' @return An object of class \code{summary.fitcox}.
#' @examples
#' summary(fitcox(remission ~ sample, delta="censor", data=leukemia))
#' @export

summary.fitcox <- function(model) {
  results <- list(coefficients=coef(model))
  class(results) <- "summary.fitcox"
  return(results)
}

#' Print a summary for a fitted Cox model
#'
#' @param summary A \code{summary.fitcox} object returned by \code{summary.fitcox}.
#' @return The object of class \code{summary.fitcox}.
#' @examples
#' print(summary(fitcox(remission ~ sample, delta="censor", data=leukemia)))
#' @export

print.summary.fitcox <- function(summary) {
  print(summary$coefficients)
  invisible(summary)
}

#' Generate predictions for a fitted Cox model
#'
#' @param model A \code{fitcox} object returned by the \code{fitcox} function.
#' @param newdata A dataframe of values to generate predictions for.
#' @param type The type of prediction to make: "median" | "expected" | "survival"
#' @return A vector of predictions for the new data
#' @details The types of prediction are:
#' \describe{
#'    \item{median}{The median survival times given the predictors}
#'    \item{expected}{The expected number of events at the event time}
#'    \item{survival}{The estimated survival probability at the event time}
#' }
#' @examples
#' predict(fitcox(remission ~ sample, delta="censor", data=leukemia), newdata=leukemia)
#' @importFrom stats predict
#' @importFrom utils tail
#' @export

predict.fitcox <- function(model, newdata, type="median") {
  formula <- model$formula

  # Add an event time of zero at the start for times before the first event
  event.times <- c(0, model$event.times)
  n <- nrow(newdata)
  results <- rep(NA, n)

  # Extract the event times for the new data
  response.var <- all.vars(formula)[1]
  responses <- newdata[[response.var]]

  # Calculate exp(z*beta) for the new data
  X <- model.matrix(formula, data=newdata)
  beta <- c(0, model$beta)
  exp.betaX <- exp(X %*% beta)

  # Baseline hazard and cumulative hazard. Cumulative hazard starts at zero.
  h0 <- model$h0
  H0 <- c(0, cumsum(h0))

  # Find the index of the event time corresponding to the response time
  closest.time.indices <- vapply(responses,
                                 function(t) tail(which(t >= event.times), n=1),
                                 numeric(1))

  if (type=="expected") {
    results <- H0[closest.time.indices]*exp.betaX
  }
  else if (type=="survival") {
    results <- exp(-H0[closest.time.indices]*exp.betaX)
  }
  else if (type=="median") {
    # Values of H(t) such that S(t) = 0.5
    H <- log(2)/exp.betaX
    # Need to find time for largest value of H0 such that H >= H0
    H.indices <- vapply(H, function(x) tail(which(x >= H0), n=1), numeric(1))
    results <- event.times[H.indices]
  }

  return(as.vector(results))
}

#' Generate concordance statistics for a fitted Cox model
#'
#' @param model A \code{fitcox} object returned by the \code{fitcox} function.
#' @param newdata A dataframe of values to generate predictions for.
#' @return A list containing:
#' \describe{
#'    \item{n.concordant}{The number of concordant pairs}
#'    \item{n.discordant}{The number of disccordant pairs}
#'    \item{c.index}{Harrel's c-index}
#' }
#' @examples
#' concordance.fitcox(fitcox(remission ~ sample, delta="censor", data=leukemia), leukemia)
#' @importFrom DescTools ConDisPairs
#' @export

concordance.fitcox <- function(model, newdata) {
  # Extract the event times for the new data
  formula <- model$formula
  response.var <- all.vars(formula)[1]
  responses <- newdata[[response.var]]

  medians <- predict(model, newdata=newdata, type="median")
  CDPairs <- ConDisPairs(table(responses, medians))
  n.C <- CDPairs$C
  n.D <- CDPairs$D
  c.index <- n.C/(n.C + n.D)
  return(list(n.concordant=n.C, n.discordant=n.D, c.index=c.index))
}

#' Generate cross-validation error for a fitted Cox model
#'
#' @param data A dataframe of values containing predictors and response.
#' @param model A \code{fitcox} object returned by the \code{fitcox} function.
#' @param K The number of groups to use for cross-validation
#' @return The mean concordance index generated using cross-validated
#' @examples
#' cv.fitcox(leukemia, fitcox(remission ~ sample, delta="censor", data=leukemia), 5)
#' @export

cv.fitcox <- function(data, model, K) {
  # Formula object and delta to use for model fitting
  formula <- model$formula
  delta <- model$delta

  cv.errors <- rep(NA, K)
  n <- nrow(data)

  # Current indices contain the indices yet to be selected
  current.indices <- 1:n

  for (i in 1:K) {
    # Generate fold indices for this value of i
    fold.size <- round(i*n/K) - round((i-1)*n/K)
    fold.indices <- sample(current.indices, fold.size)
    current.indices <- current.indices[! current.indices %in% fold.indices]

    # Generate train and test data sets
    train.data <- data[-fold.indices, ]
    test.data <- data[fold.indices, ]

    # Train model on the training data, and concordance index on test data
    train.model <- fitcox(formula, delta, train.data)
    cv.errors[i] <- concordance.fitcox(train.model, test.data)$c.index
  }

  return(mean(cv.errors))
}
