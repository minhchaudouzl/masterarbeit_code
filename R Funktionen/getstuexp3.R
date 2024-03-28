library(BBmisc)

# Stückweise Exponentialfunktion mit 3 Parametern

getstuexp3 <- function (p = c(0.025, 0.5, 0.975), q, show.output = TRUE, plot = TRUE, 
                        start,
                        tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 
                                                                                 0.9), wert1, wert2,...) 
{
  if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
    stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", 
         call. = FALSE)
  }
  if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == 
                                                     seq(1:length(q))) == 0) {
    stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", 
         call. = FALSE)
  }
  if (min(p) < 0 | max(p) > 1) {
    stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", 
         call. = FALSE)
  }
  if (min(q) < 0) {
    stop("INVALID INPUT, percentiles are out of the domain [0, inf) => Weibull distribution couldn't be fitted!", 
         call. = FALSE)
  }
  if (length(p) != length(q) | length(p) != length(fit.weights) | 
      length(q) != length(fit.weights)) {
    stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", 
         call. = FALSE)
  }
  if (length(q) < 2) {
    stop("INVALID INPUT, at least two quantiles must be known!", 
         call. = FALSE)
  }
  if (!is.logical(show.output)) {
    stop("INVALID INPUT, the argument 'show.output' should be logical!", 
         call. = FALSE)
  }
  if (!is.logical(plot)) {
    stop("INVALID INPUT, the argument 'plot' should be logical!", 
         call. = FALSE)
  }
  if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
    stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", 
         call. = FALSE)
  }
  fit.weights.original <- fit.weights
  fit.weights <- fit.weights/sum(fit.weights)
  minimize <- function(theta) {
    summand <- suppressWarnings(((q > 0 & q <= wert1 ) * (1 - exp(-theta[1] * q)) + 
                                   (q > wert1 & q <= wert2 ) *  (1 - exp(-theta[1] * wert1))  + 
                                   (exp(-wert1 * theta[2]) - exp(-theta[2] * q)) +
                                   (q > wert2 & q <= 65 ) * (1 - exp(-theta[1] * wert1)) + 
                                   (exp(-wert1 * theta[2]) - exp(-wert2 * theta[2])) + 
                                   (exp(-wert2 * theta[3]) - exp(-theta[3] * q))) - p)
                                  
    summand <- summand * fit.weights
    sum(summand^2)
  }
  fit <- c()
  fit$value <- tol + 1
  if(is.null(start)) start <- c(0,1)
  try1 <- try(fit <- stats::optim(par = start, minimize, 
                                  method = "L-BFGS-B", lower = c(0.001, 0.001), upper = c(10000, 
                                                                                          10000)), silent = TRUE)
  if (is.error(try1) || fit$value >= tol) {
    warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", 
            call. = FALSE)
    fit <- c()
    fit$value <- tol + 1
    try2 <- try(fit <- stats::optim(par = start, minimize, 
                                    method = "BFGS"), silent = TRUE)
    if (is.error(try2) || fit$value >= tol) {
      warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", 
              call. = FALSE)
      Par <- NA
    }
    else if (fit$value < tol) {
      message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)")
      Par <- fit$par
      names(Par) <- c("1.para", "2.para", "3.para")
      if (show.output) 
        print(fit)
    }
  }
  else if (fit$value < tol) {
    message("The fitting procedure 'L-BFGS-B' was successful!")
    Par <- fit$par
    names(Par) <- c("1.para", "2.para", "3.para")
    if (show.output) 
      print(fit)
  }
  if (prod(!is.na(Par)) & plot) {
    main1 <- paste("1.para = ", signif(Par["1.para"], 
                                     digits = 3))
    main2 <- paste("2.para = ", signif(Par["2.para"], 
                                       digits = 3)) 
    main3 <- paste("3.para = ", signif(Par["3.para"], 
                                       digits = 3))
    main <- paste("3 stueckw. Exponential (", main1, ", ", main2, ", ", main3, ")", sep = "")
    
    sub = paste("fit.weights = c(", paste(fit.weights.original, 
                                          collapse = ", "), ")", sep = "")
    Support.lim <- c(0, 60)
    Support <- seq(Support.lim[1], Support.lim[2], length = 200)
    Probability <- (((Support > 0 & Support <= wert1 ) * (1 - exp(-Par["1.para"] * Support)) + 
                       (Support > wert1 & Support <= wert2 ) * (1 - exp(-Par["1.para"] * wert1)) +
                                                          (exp(-wert1 * Par["2.para"]) - exp(-Par["2.para"] * Support)) +
                       (Support > wert2 & Support <= 65 ) * (1 - exp(-Par["1.para"] * wert1)) + 
                       (exp(-wert1 * Par["2.para"]) - exp(-wert2 * Par["2.para"])) + 
                       (exp(-wert2 * Par["3.para"]) - exp(-Par["3.para"] * Support))))
                                                      
                                                      
    graphics::plot(Support, Probability, type = "l", 
                   xlim = range(Support.lim, q), main = main, xlab = "Quantiles", 
                   sub = sub, ...)
    graphics::points(x = q, y = p, pch = 19, ...)
  }
  return(Par)
}