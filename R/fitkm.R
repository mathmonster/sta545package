#' Kaplan-Meier product limit estimator
#'
#' \code{fitkm} fits a Kaplan-Meier product limit estimator.
#'
#' @param formula A formula object with the time variable on the left-hand side.
#' @param delta Censoring indicator, 0 = right-censored, 1 = event occurred.
#' @param data A data frame containing the variables referred to in \code{formula}.
#' @return An object of class \code{fitkm}, containing one or more
#'     \code{\link{fitkm.curve}} objects.
#' @examples
#' fitkm(recovery ~ 1, delta="status", data=data1)
#' fitkm(remission ~ sample, delta="censor", data=leukemia)
#' @seealso \code{\link{summary.fitkm}}, \code{\link{plot.fitkm}},
#'     \code{\link{fitkm.curve}}
#' @export

fitkm <- function(formula, delta, data) {
  # Extract the response and predictor variable (if a predictor is used)
  response.var <- all.vars(formula)[1]
  predictor <- all.vars(formula)[2]
  # NA predictor means use all the data to construct a single survival curve
  # Otherwise, partition the data by levels of the predictor
  if (is.na(predictor)) {
    partitions <- list(alldata=data)
  }
  else {
    partitions <- split(data, data[[predictor]])
  }

  results <- list()

  # Generate a separate KM curve for each level of the predictor variable
  for (level in names(partitions)) {
    subdata <- partitions[[level]]
    n <- nrow(subdata)
    deltas <- subdata[[delta]]
    all.times <- subdata[[response.var]]
    events <- all.times[deltas==1]
    event.times <- sort(unique(events))
    censor.times <- sort(unique(all.times[deltas==0]))

    # Calculate the number of failures at each event time
    n.events <- vapply(event.times, function(x) sum(x==events), numeric(1))

    # Calculate the size of the risk set at each event time
    n.risk <- vapply(event.times, function(x) sum(x <= all.times), numeric(1))

    # Calculate the survival function at each event time
    surv.prob <- cumprod(1 - n.events/n.risk)

    # Return results in an object of class "fitkm.object"
    curve <- list(n=n, level=level, time=event.times, n.events=n.events,
                  n.risk=n.risk, surv.prob=surv.prob, censored=censor.times)
    class(curve) <- "fitkm.curve"
    results[[level]] <- curve
  }

  class(results) <- "fitkm"
  return(results)
}

#' Print details of a set of one or more Kaplan-Meier curves
#'
#' @param model A \code{fitkm} object returned by the \code{\link{fitkm}} function.
#' @return The object of class \code{fitkm}.
#' @examples
#' print(fitkm(recovery ~ 1, delta="status", data=data1))
#' @export

print.fitkm <- function(model) {
  print.data <- data.frame(n=vapply(model, function(x) x$n, numeric(1)),
                           n.events=vapply(model, function(x) sum(x$n.events), numeric(1))
                           )
  print(print.data)
  invisible(model)
}

#' Summarize a set of one or more Kaplan-Meier curves
#'
#' @param model A \code{fitkm} object returned by the \code{\link{fitkm}} function.
#' @return A list of summaries for each \code{fitkm.curve} object
#' @examples
#' summary(fitkm(recovery ~ 1, delta="status", data=data1))
#' @export

summary.fitkm <- function(model) {
  return(lapply(model, summary))
}

#' Produces a plot for a set of one or more Kaplan-Meier curves
#'
#' @param model A \code{fitkm} object returned by the \code{\link{fitkm}} function.
#' @param xlab String for the x-axis label
#' @param ylab String for the y-axis label
#' @param legend.cex Font size for the legend
#' @param ... Further graphical parameters to be passed to \code{plot}
#' @return An object of class \code{summary.fitkm}.
#' @examples
#' plot(fitkm(recovery ~ 1, delta="status", data=data1))
#' plot(fitkm(remission ~ sample, delta="censor", data=leukemia))
#' @importFrom graphics lines legend
#' @importFrom grDevices hcl.colors
#' @export

plot.fitkm <- function(model, xlab="Time", ylab="Survival Probability", legend.cex=1, ...) {
  n.curves <- length(model)
  if (n.curves > 1) {
    line.colors <- hcl.colors(n.curves, palette="Zissou 1")
  }
  else {
    line.colors <- c("red")
  }

  # Plot the first curve
  plot(model[[1]], xlab=xlab, ylab=ylab, col=line.colors[1],  ylim=c(0, 1), ...)

  # Add lines for the remaining curves
  if (n.curves > 1) {
    for (i in 1:n.curves) {
      curve <- model[[i]]
      lines(curve, col=line.colors[i], ...)
      }
    }

  # Add a legend if more than one curve being plotted
  if (n.curves > 1) {
    legend("topright", legend=names(model), col=line.colors, lwd=2, cex=legend.cex)
  }
}

#' Fitted Kaplan-Meier survival curve object
#'
#' A survival curve object containing the results of fitting a Kaplan-Meier
#' survival curve. The notation corresponds to that on page 188 of the Cox 1972
#' paper
#'
#' @param n the number of observations
#' @param level the factor level that the survival curve has been generated for
#' @param time t(i), the unique failure times sorted in ascending order
#' @param n.events  m(i) , the number of failure times at each t(i)
#' @param n.risk r(i), the number in the risk set at each t(i)
#' @param surv.prob F(i), the estimated survival function at each t(i)
#' @param censored the unique censoring times sorted in ascending order
#'
#' @seealso \code{\link{fitkm}}
#' @export

fitkm.curve <- function(n, level, time, n.events, n.risk, surv.prob, censored) {
}

#' Print details of a fitted Kaplan-Meier curve
#'
#' @param curve A \code{fitkm.curve} object returned by the \code{fitkm} function.
#' @return An object of class \code{summary.fitkm}.
#' @export

print.fitkm.curve <- function(curve) {
  print(data.frame(level=curve$level, n=curve$n, n.events=sum(curve$n.events)))
  invisible(curve)
}

#' Produces a summary for a fitted Kaplan-Meier curve
#'
#' @param curve A \code{fitkm.curve} object returned by the \code{\link{fitkm}} function.
#' @return A dataframe containing summary information for each event time.
#' @examples
#' summary(fitkm(recovery ~ 1, delta="status", data=data1))
#' @export

summary.fitkm.curve <- function(curve) {
  return(data.frame(time=curve$time,
                    n.events=curve$n.events,
                    n.risk=curve$n.risk,
                    surv.prob=curve$surv.prob))
}

#' Produces a plot for a single fitted Kaplan-Meier curve
#'
#' @param curve A \code{fitkm.curve} object returned by the \code{fitkm} function.
#' @param censor Plot the censored observations if \code{TRUE}
#' @param ... Further graphical parameters
#' @importFrom graphics points
#' @export

plot.fitkm.curve <- function(curve,  censor=TRUE, ...) {
  times <- c(0, curve$time)
  surv.prob <- c(1, curve$surv.prob)
  plot(times, surv.prob, type="s",  ...)

  # Plot the censored observations
  if (censor) {
    censored <- curve$censored
    censored.s <- vapply(censored, function(t) surv.prob[max(which(t >= times))], numeric(1))
    points(censored, censored.s, type="p", pch="|", ...)
  }
}

#' Plots a line for a single fitted Kaplan-Meier curve
#'
#' @param curve A \code{fitkm.curve} object returned by the \code{fitkm} function.
#' @param censor Plot the censored observations if \code{TRUE}
#' @param ... Further graphical parameters
#' @importFrom graphics lines points
#' @export

lines.fitkm.curve <- function(curve, censor=TRUE,  ...) {
  times <- c(0, curve$time)
  surv.prob <- c(1, curve$surv.prob)
  lines(times, surv.prob, type="s",  ...)

  # Plot the censored observations
  if (censor) {
    censored <- curve$censored
    censored.s <- vapply(censored, function(t) surv.prob[max(which(t >= times))], numeric(1))
    points(censored, censored.s, type="p", pch="|", ...)

  }
}
