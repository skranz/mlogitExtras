#' Draw random coefficients of individuals for an estimated mixed logit model
#'
#' Note that the function is still restricted and does not yet work for all sorts
#' of random parameter distributions that mlogit supports.
#'
#' @param mod the estimated models
#' @param num.draws the number of individuals for which random coefficients shall be drawn
#' @param use.halton shall Halton sequences be used instead of psedudo-random numbers?
#' @return A matrix with num.draws rows and one column for each explanatory variable that has random coefficients
ml_draw_betas = function(mod, num.draws=100, use.halton = FALSE) {
  restore.point("ml_draw_betas")

  rpar=mod$rpar
  rvars = names(rpar)

  # Draw individual betas
  beta.mat = matrix(NA, nrow=num.draws, ncol=length(rvars))
  colnames(beta.mat) = rvars

  if (use.halton) {
    library(randtoolbox)
    # Use halton draws for beta.mat
    beta.unif = randtoolbox::halton(num.draws, dim=length(rvars))
    for (k in seq_along(rvars)) {
      beta.mat[, k] = transform.uniform(beta.unif[,k],mod$rpar[[rvars[k]]])
    }
  } else {
    # Use pseudo-random numbers for beta.mat
    for (k in seq_along(rvars)) {
      #beta.mat[, k] = draw.iid.rpar.beta(mod$rpar[[rvars[k]]], num.draws=num.draws)
      beta.mat[, k] = transform.uniform(runif(num.draws),mod$rpar[[rvars[k]]])
    }
  }
  beta.mat
}

draw.iid.rpar.beta = function(rpar, num.draws=1000) {
  if (rpar$dist == "n") {
    rnorm(num.draws,rpar$mean, abs(rpar$sigma))
  } else if (rpar$dist == "u") {
    half.span = rpar$sigma / 2
    runif(num.draws, rpar$mean-half.span, rpar$mean+half.span)
  } else {
    stop("Prediction is implemented so far only for iid normal and iid uniform distributed random coefficients.")
  }
}

# Transform uniformely distributed variables
transform.uniform = function(beta.unif, rpar) {
  if (rpar$dist == "n") {
    qnorm(beta.unif, mean = rpar$mean, sd = abs(rpar$sigma))
  } else if (rpar$dist == "u") {
    span = rpar$sigma
    beta.unif*span+(rpar$mean-span/2)
  } else if (rpar$dist == "t") {
    min = rpar$mean-0.5*abs(rpar$sigma)
    max = rpar$mean+0.5*abs(rpar$sigma)
    qtri(beta.unif, min,max)
  } else if (rpar$dist == "ln") {
    qlnorm(beta.unif, mean = rpar$mean, sd = abs(rpar$sigma))
  } else if (rpar$dist == "cn") {
    beta = qnorm(beta.unif, mean = rpar$mean, sd = abs(rpar$sigma))
    beta[beta<0] = 0
    beta
  } else {
    stop("Prediction is implemented so far only for iid normal and iid uniform distributed random coefficients.")
  }


}

